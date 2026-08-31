#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-04-19 18:57

import argparse
import importlib.metadata
import json
import logging
import sys
from dataclasses import dataclass
from typing import Optional

import cutadapt
from cutadapt.files import InputPaths, OutputFiles
from cutadapt.pipeline import PairedEndPipeline, SingleEndPipeline
from cutadapt.predicates import Predicate, TooManyN, TooShort
from cutadapt.report import Statistics, minimal_report
from cutadapt.runners import make_runner
from cutadapt.steps import (
    PairedEndFilter,
    PairedEndSink,
    SingleEndFilter,
    SingleEndSink,
)
from cutadapt.modifiers import (
    AdapterCutter,
    PairedEndRenamer,
    QualityTrimmer,
    Renamer,
    SuffixRemover,
    SingleEndModifier,
)
from cutadapt.utils import Progress

from .common import (
    BUILDIN_ADAPTERS,
    print_builtin_adapters,
    remove_fq_suffix,
)
from .grammar import (
    CompiledScheme,
    load_scheme_file,
    parse_scheme,
)


class ReverseComplementModifier(SingleEndModifier):
    """Reverse-complement a read (used with ``--auto-rc`` for ``-`` libraries)."""

    def __call__(self, read, info):
        return read.reverse_complement()


class IsUntrimmedAny(Predicate):
    """Select reads where at least one expected inline-barcode adapter was not
    found in ``info.matches`` — i.e. reads that missed an inline barcode.
    Used by ``--ensure-inline-barcode`` to route such reads to the discard
    output. ``ref_adapters`` are the exact Adapter objects compiled into the
    5' inline-barcode steps, compared by identity.
    """

    def __init__(self, ref_adapters):
        self.ref_adapters = list(ref_adapters)

    def __repr__(self):
        return f"IsUntrimmedAny(ref_adapters={self.ref_adapters!r})"

    def test(self, read, info):
        matched = [m.adapter for m in info.matches]
        return any(a not in matched for a in self.ref_adapters)


class LowAverageQuality(Predicate):
    """Select reads whose mean Phred quality is below a threshold."""

    def __init__(self, min_avg_quality):
        self.min_avg_quality = min_avg_quality

    def __repr__(self):
        return f"LowAverageQuality(min={self.min_avg_quality})"

    def test(self, read, info):
        if len(read) == 0 or read.qualities is None:
            return False
        # Phred+33 encoding: each byte is ord(qual) - 33.
        return (sum(b - 33 for b in read.qualities.encode()) / len(read)
                < self.min_avg_quality)


class _TaggedSingleEndFilter(SingleEndFilter):
    """SingleEndFilter that appends ``reason=...`` to the read name before
    writing, so a single shared discard file can carry multiple reasons."""

    def __init__(self, predicate, writer, reason):
        super().__init__(predicate, writer)
        self._reason = reason

    def __call__(self, read, info):
        if self._predicate.test(read, info):
            self._filtered += 1
            if self._writer is not None:
                read.name = f"{read.name} reason={self._reason}"
                self._writer.write(read)
            return None
        return read


class _TaggedPairedEndFilter(PairedEndFilter):
    """PairedEndFilter that appends ``reason=...`` to both read names before
    writing, so a single shared discard file can carry multiple reasons."""

    def __init__(self, predicate1, predicate2, writer, reason,
                 pair_filter_mode="any"):
        super().__init__(predicate1, predicate2, writer, pair_filter_mode)
        self._reason = reason

    def __call__(self, read1, read2, info1, info2):
        if self._is_filtered(read1, read2, info1, info2):
            self._filtered += 1
            if self.writer is not None:
                read1.name = f"{read1.name} reason={self._reason}"
                read2.name = f"{read2.name} reason={self._reason}"
                self.writer.write(read1, read2)
            return None
        return read1, read2

#  statistics workaround ....
# cutadapt's Statistics._collect_modifier contains an assert that can be
# violated by exotic modifier orderings. Guard it at run time (in main()),
# never at import time, and log instead of silently swallowing.
_original_collect_modifier = Statistics._collect_modifier


def _guarded_collect_modifier(self, *args, **kwargs):
    try:
        return _original_collect_modifier(self, *args, **kwargs)
    except AssertionError:
        logging.warning(
            "cutadapt Statistics._collect_modifier raised AssertionError; "
            "statistics may be incomplete. This is a cutadapt quirk, not a "
            "cutseq bug."
        )
        return None


def _apply_statistics_guard():
    """Scope the Statistics guard to run time (see module docstring)."""
    if Statistics._collect_modifier is _original_collect_modifier:
        Statistics._collect_modifier = _guarded_collect_modifier


def _option_arities(parser):
    """Map each option string to how many values it consumes.

    Derived directly from the parser's actions so it cannot drift from the CLI
    spec: ``0`` for flags (``store_true`` / ``BooleanOptionalAction``), ``1``
    for single-value options, and ``'many'`` for the varargs ``-o``/``-d``
    options.
    """
    arities = {}
    for opt, action in parser._option_string_actions.items():
        n = action.nargs  # None -> 1; 0 -> 0; "+" -> many
        arities[opt] = "many" if n == "+" else (1 if n is None else 0)
    return arities


def _reorder_inputs_first(argv, arities):
    """Reorder ``argv`` so the positional input reads come before the options.

    argparse greedily assigns *every* trailing token to an ``nargs='+'`` option
    (``-o`` / ``-d``) when the ``nargs='*'`` ``input_file`` positional is still
    empty, so ``cutseq -A X -o a.fq b.fq in.fq`` swallows ``in.fq`` as an
    output file and then fails with "Input file is required". Reordering lets
    input files be given before or after the options.

    The two valid input counts (2 = paired, then 1 = single-end) are tried as
    a fixed point: a count ``k`` is consistent when the walk consumes exactly
    ``k`` values for each varargs option and collects exactly ``k`` leftover
    positional tokens. The first consistent count wins; argv is rewritten with
    those positionals moved to the front. If no count is consistent (a genuine
    user error, or ``--`` is used so argparse should handle it), argv is
    returned unchanged and argparse reports its own error.
    """
    if not argv or "--" in argv:
        return argv
    # If a real positional token appears before the first varargs option
    # (-o/-d), argparse can already disambiguate; leave argv untouched.
    # Option *values* (e.g. 'SMALLRNA' after -A) are skipped during the walk.
    i, n = 0, len(argv)
    identity_safe = False
    while i < n:
        t = argv[i]
        if t in arities:
            a = arities[t]
            if a == "many":
                break
            i += 1 + (1 if a == 1 else 0)
        else:
            identity_safe = True
            break
    if identity_safe:
        return argv
    for k in (2, 1, 0):
        inputs, rest = [], []
        i, n = 0, len(argv)
        ok = True
        while i < n and ok:
            tok = argv[i]
            if tok in arities:
                arity = arities[tok]
                rest.append(tok)
                i += 1
                if arity == "many":
                    for _ in range(k):
                        if i >= n or argv[i] in arities:
                            ok = False
                            break
                        rest.append(argv[i])
                        i += 1
                elif arity == 1:
                    if i >= n:
                        ok = False
                        break
                    rest.append(argv[i])
                    i += 1
                # arity == 0: a flag, nothing to consume
            else:
                inputs.append(tok)
                i += 1
        if ok and len(inputs) == k:
            return inputs + rest
    return argv


__version__ = importlib.metadata.version(__package__ or __name__)


@dataclass
class CutadaptConfig:
    """Centralized settings for the cutadapt pipeline.

    Values mirror the CLI defaults so the pipelines can be driven directly
    (tests, library use) without argparse.
    """

    rname_suffix: bool = False
    ensure_inline_barcode: bool = False
    conditional_cutter: bool = True
    min_length: int = 20
    min_quality: int = 20
    min_avg_quality: Optional[float] = None
    max_n: Optional[float] = None
    auto_rc: bool = False
    dry_run: bool = False
    threads: int = 1
    json_file: Optional[str] = None
    force_trim_min_length: int = 50
    force_anywhere: bool = False
    name_format: Optional[str] = None
    auto_inline: bool = True


def json_report(
    file,
    stats,
    barcode,
    input1,
    input2,
    output1,
    output2,
    discard1,
    discard2,
):
    """
    Generates a JSON report summarizing the trimming statistics.

    The report includes cutadapt version, input/output file paths, barcode configuration,
    and detailed trimming statistics from the cutadapt run.

    :param file: Path to the output JSON file.
    :type file: str
    :param stats: Cutadapt Statistics object.
    :type stats: cutadapt.report.Statistics
    :param barcode: Dict describing the adapter scheme captures used for the run.
    :type barcode: dict
    :param input1: Path to the first input FASTQ file.
    :type input1: str
    :param input2: Path to the second input FASTQ file (None for single-end).
    :type input2: str, optional
    :param output1: Path to the first output trimmed FASTQ file.
    :type output1: str
    :param output2: Path to the second output trimmed FASTQ file (None for single-end).
    :type output2: str, optional
    :param discard1: Path to the first discard output (reason-tagged reads; None if disabled).
    :type discard1: str, optional
    :param discard2: Path to the second discard output (None if disabled).
    :type discard2: str, optional
    """
    d = {
        "tag": "Cutadapt report",
        "cutadapt_version": cutadapt.__version__,
        "input": {
            "path1": input1,
            "path2": input2,
            "paired": True
            if input2
            else False,  # Corrected based on presence of input2
        },
        "output": {
            "output1": output1,
            "output2": output2,
            "discard1": discard1,
            "discard2": discard2,
        },
        "barcode": barcode,
    }

    d.update(stats.as_json())
    # adapters_read1/2 -> trimmed_lengths to emppty list
    # This is done to reduce the size of the JSON report, as detailed trimmed_lengths
    # can be very large and are often not needed for summary purposes.
    if d.get("adapters_read1"):
        for m in d["adapters_read1"]:
            if m.get("five_prime_end"):
                m["five_prime_end"]["trimmed_lengths"] = []
            if m.get("three_prime_end"):
                m["three_prime_end"]["trimmed_lengths"] = []
    if d.get("adapters_read2"):
        for m in d["adapters_read2"]:
            if m.get("five_prime_end"):
                m["five_prime_end"]["trimmed_lengths"] = []
            if m.get("three_prime_end"):
                m["three_prime_end"]["trimmed_lengths"] = []
    with open(file, "w") as json_file:
        json_file.write(json.dumps(d, indent=2))


def _resolve_scheme(scheme, auto_inline=True):
    """Return ``(orientation, left, right)`` for a scheme value (string/YAML/tokens)."""
    if isinstance(scheme, tuple):
        return scheme
    is_yaml, resolved = load_scheme_file(scheme)
    if is_yaml:
        return resolved
    return parse_scheme(resolved, auto_inline=auto_inline)


_renamers = (Renamer, PairedEndRenamer)


def _print_known_primers():
    """Print the known sequencing primers used for inline-barcode detection."""
    from .primers import MIN_PRIMER_MATCH, SEQUENCING_PRIMERS

    print(f"\nKnown sequencing primers (terminal match >= {MIN_PRIMER_MATCH} bp, "
          f"either strand) used for inline-barcode auto-detection:\n")
    width = max(len(n) for n in SEQUENCING_PRIMERS)
    for name, seq in SEQUENCING_PRIMERS.items():
        print(f"{name.ljust(width)}  {seq}")
    print("\nUse --no-auto-inline to disable auto-detection.\n")


def _describe(mod):
    """Render a modifier readably for dry-run logs (native reprs are opaque)."""
    if isinstance(mod, tuple):
        return f"({_describe(mod[0])}, {_describe(mod[1])})"
    if isinstance(mod, AdapterCutter):
        parts = ", ".join(a.sequence for a in mod.adapters)
        return f"AdapterCutter({parts})"
    if isinstance(mod, _renamers):
        return f"{type(mod).__name__}({mod._template!r})"
    if hasattr(mod, "_template"):
        name = getattr(mod, "__name__", type(mod).__name__)
        if name == "function":  # fast-renamer closure
            name = "Renamer"
        return f"{name}({mod._template!r})"
    return repr(mod)


def _build_scheme(scheme, settings):
    """Resolve a scheme value (built-in name / grammar string / YAML tokens)
    into a CompiledScheme bound to the current pipeline settings."""
    orientation, left, right = _resolve_scheme(scheme,
                                               auto_inline=settings.auto_inline)
    return CompiledScheme(
        orientation, left, right,
        conditional_cutter=settings.conditional_cutter,
        force_trim_min_length=settings.force_trim_min_length,
        force_anywhere=settings.force_anywhere,
    )


def _scheme_modifiers(cs, paired, settings):
    """Build the shared modifier list (trimming steps) for a CompiledScheme."""
    if paired:
        mods1, mods2 = cs.modifiers(paired=True)
        n = max(len(mods1), len(mods2))
        modifiers = []
        if settings.rname_suffix:
            modifiers.append((SuffixRemover(".1"), SuffixRemover(".2")))
            modifiers.append((SuffixRemover("/1"), SuffixRemover("/2")))
        for i in range(n):
            modifiers.append((
                mods1[i] if i < len(mods1) else None,
                mods2[i] if i < len(mods2) else None,
            ))
        modifiers.append((
            QualityTrimmer(cutoff_front=0, cutoff_back=settings.min_quality),
            QualityTrimmer(cutoff_front=0, cutoff_back=settings.min_quality),
        ))
        if settings.auto_rc:
            logging.warning(
                "--auto-rc is ignored for paired-end data: R1/R2 are already "
                "oriented via the scheme's strand marker."
            )
        modifiers.append(cs.renamer(paired=True,
                                    name_format=settings.name_format))
        return modifiers

    modifiers, _ = cs.modifiers(paired=False)
    if settings.rname_suffix:
        modifiers = [SuffixRemover(".1"), SuffixRemover("/1")] + modifiers
    modifiers = modifiers + [QualityTrimmer(cutoff_front=0,
                                            cutoff_back=settings.min_quality)]
    if settings.auto_rc:
        if cs.orientation == "-":
            modifiers.append(ReverseComplementModifier())
        else:
            logging.warning(
                "Library scheme is not '-' strand, but --auto-rc is enabled; "
                "ignoring --auto-rc."
            )
    modifiers = modifiers + [cs.renamer(paired=False,
                                        name_format=settings.name_format)]
    return modifiers


def pipeline_grammar_single(input1, output1, discard1, scheme, settings):
    """Run a single-end library-grammar pipeline (see cutseq.grammar).

    ``discard1`` is the path for discarded reads; each discarded read has its
    name tagged with the reason (``reason=too_short``, ``reason=too_many_n``,
    ``reason=low_quality`` or ``reason=no_barcode``). Pass None to disable the
    discard output.
    """
    cs = _build_scheme(scheme, settings)
    modifiers = _scheme_modifiers(cs, paired=False, settings=settings)

    if settings.dry_run:
        for i, m in enumerate(modifiers, 1):
            print(f"Step {i}: {_describe(m)}")
        return

    inpaths = InputPaths(input1)
    with make_runner(inpaths, cores=settings.threads) as runner:
        outfiles = OutputFiles(
            proxied=settings.threads > 1,
            qualities=runner.input_file_format().has_qualities(),
            interleaved=False,
        )
        steps = []
        discard_writer = (outfiles.open_record_writer(discard1, interleaved=False)
                          if discard1 is not None else None)
        if settings.max_n is not None:
            steps.append(_TaggedSingleEndFilter(
                TooManyN(settings.max_n), discard_writer, "too_many_n"))
        if settings.min_avg_quality is not None:
            steps.append(_TaggedSingleEndFilter(
                LowAverageQuality(settings.min_avg_quality), discard_writer,
                "low_quality"))
        steps.append(_TaggedSingleEndFilter(
            TooShort(settings.min_length), discard_writer, "too_short"))
        r1_adps, _ = cs.inline_adapters(paired=False)
        if settings.ensure_inline_barcode and r1_adps and discard_writer is not None:
            steps.append(_TaggedSingleEndFilter(
                IsUntrimmedAny(r1_adps), discard_writer, "no_barcode"))
        steps.append(SingleEndSink(
            outfiles.open_record_writer(output1, interleaved=False),
        ))
        pipeline = SingleEndPipeline(modifiers, steps)
        stats = runner.run(pipeline, Progress(), outfiles)

        if settings.json_file is not None:
            json_report(settings.json_file, stats, cs.summary(), input1, None,
                        output1, None, discard1, None)
        print(minimal_report(stats, time=None, gc_content=None), file=sys.stderr)
    outfiles.close()


def pipeline_grammar_paired(input1, input2, output1, output2, discard1, discard2,
                            scheme, settings):
    """Run a paired-end library-grammar pipeline (see cutseq.grammar).

    ``discard1``/``discard2`` are the paths for discarded read pairs; each
    discarded pair has both read names tagged with the reason (``reason=too_short``,
    ``reason=too_many_n``, ``reason=low_quality`` or ``reason=no_barcode``).
    Pass None to disable the discard output.
    """
    cs = _build_scheme(scheme, settings)
    modifiers = _scheme_modifiers(cs, paired=True, settings=settings)

    if settings.dry_run:
        for i, b in enumerate(modifiers, 1):
            print(f"Step {i}: {_describe(b)}")
        return

    inpaths = InputPaths(input1, input2)
    with make_runner(inpaths, cores=settings.threads) as runner:
        outfiles = OutputFiles(
            proxied=settings.threads > 1,
            qualities=runner.input_file_format().has_qualities(),
            interleaved=False,
        )
        steps = []
        discard_writer = (outfiles.open_record_writer(discard1, discard2,
                                                      interleaved=False)
                          if (discard1 is not None and discard2 is not None)
                          else None)
        if settings.max_n is not None:
            steps.append(_TaggedPairedEndFilter(
                TooManyN(settings.max_n), TooManyN(settings.max_n),
                discard_writer, "too_many_n"))
        if settings.min_avg_quality is not None:
            steps.append(_TaggedPairedEndFilter(
                LowAverageQuality(settings.min_avg_quality),
                LowAverageQuality(settings.min_avg_quality),
                discard_writer, "low_quality"))
        steps.append(_TaggedPairedEndFilter(
            TooShort(settings.min_length), TooShort(settings.min_length),
            discard_writer, "too_short"))
        r1_adps, r2_adps = cs.inline_adapters(paired=True)
        if (settings.ensure_inline_barcode and discard_writer is not None
                and (r1_adps or r2_adps)):
            steps.append(_TaggedPairedEndFilter(
                IsUntrimmedAny(r1_adps) if r1_adps else None,
                IsUntrimmedAny(r2_adps) if r2_adps else None,
                discard_writer, "no_barcode", pair_filter_mode="any"))
        steps.append(PairedEndSink(
            outfiles.open_record_writer(output1, output2),
        ))
        pipeline = PairedEndPipeline(modifiers, steps)
        stats = runner.run(pipeline, Progress(), outfiles)

        if settings.json_file is not None:
            json_report(settings.json_file, stats, cs.summary(), input1, input2,
                        output1, output2, discard1, discard2)
        print(minimal_report(stats, time=None, gc_content=None), file=sys.stderr)
    outfiles.close()


def run_cutseq(args):
    """
    Sets up configurations and executes the grammar-based cutadapt pipeline.

    This function initializes the CutadaptConfig based on the command-line
    arguments. It then determines whether to run the single-end or
    paired-end grammar pipeline based on the number of input files provided.

    :param args: Parsed command-line arguments from argparse.
    :type args: argparse.Namespace
    """
    settings = CutadaptConfig()
    settings.rname_suffix = args.with_rname_suffix
    settings.ensure_inline_barcode = args.ensure_inline_barcode
    settings.conditional_cutter = args.conditional_cutter
    settings.threads = args.threads
    settings.min_length = args.min_length
    settings.min_quality = args.min_quality
    settings.min_avg_quality = args.min_avg_quality
    settings.max_n = args.max_n
    settings.dry_run = args.dry_run
    settings.auto_rc = args.auto_rc
    settings.json_file = args.json_file
    settings.force_trim_min_length = args.force_trim_min_length  # Pass from args
    settings.force_anywhere = args.force_anywhere
    settings.name_format = args.name_format
    settings.auto_inline = not args.no_auto_inline
    if len(args.input_file) == 1:
        pipeline_grammar_single(
            args.input_file[0], args.output_file[0], args.discard_file[0],
            args.adapter_scheme, settings,
        )
    else:
        pipeline_grammar_paired(
            args.input_file[0], args.input_file[1],
            args.output_file[0], args.output_file[1],
            args.discard_file[0], args.discard_file[1],
            args.adapter_scheme, settings,
        )


def main():
    """
    Parses command-line arguments and initiates the adapter trimming process.

    This is the main entry point for the `cutseq run` script. It sets up
    the argument parser, validates arguments, prepares output file names,
    and then calls `run_cutseq` to perform the trimming.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s -  %(levelname)s - %(message)s",
    )
    _apply_statistics_guard()
    parser = argparse.ArgumentParser(
        description="Trim sequencing adapters from NGS data automatically using cutadapt's Python API."
    )
    # input file can be one or two for single or paired-end reads, but can not be more than two
    parser.add_argument(
        "input_file",
        type=str,
        nargs="*",
        help="Input file path for NGS data, one or two files (for single or paired-end reads).",
    )
    # output file can be number of files matching the input files, if not provided it will generate based on the output prefix,
    # if no output prefix provided it will generate based on the input file name
    parser.add_argument(
        "-A",
        "-a",
        "--adapter-scheme",
        dest="adapter_scheme",
        type=str,
        help="Adapter scheme. Either a built-in adapter name (e.g. TAKARAV3) or a "
        "custom grammar scheme string. The scheme is a top-strand molecular map "
        "written 5' to 3' (left of `+`/`-`/`:` applies to R1 as-is; right of it "
        "is the 3' continuation, matched on R2 as its reverse complement). "
        "Uppercase ACGT.. are adapters, lowercase "
        "acgt.. are inline barcodes, N.. are UMI captures, X.. are masks, and "
        "`+`, `-`, `:` split the library into R1 | R2 (sense / antisense / "
        "unstranded).",
    )
    parser.add_argument(
        "-O",
        "--output-prefix",
        type=str,
        help="Output file prefix for trimmed and discard data. "
        "If not provided, output filenames are derived from input filenames.",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        nargs="+",
        help="Output file path(s) for successfully trimmed reads. Must match number of input files.",
    )
    parser.add_argument(
        "-d",
        "--discard-file",
        type=str,
        nargs="+",
        help="Output file path(s) for discarded reads. Discarded reads carry a "
        "'reason=...' tag in their read name: 'too_short' (shorter than "
        "--min-length after trimming), 'too_many_n' (exceeds --max-n), "
        "'low_quality' (mean Phred quality below --min-avg-quality), or "
        "'no_barcode' (missing an expected inline barcode, only with "
        "--ensure-inline-barcode). Must match number of input files.",
    )
    parser.add_argument(
        "--json-file",
        type=str,
        help="Output JSON file for trimming statistics.",
    )

    parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=20,
        help="Minimum quality score for trimming read tails. (Default: 20)",
    )
    parser.add_argument(
        "-m",
        "--min-length",
        type=int,
        default=20,
        help="Minimum length of reads to keep after trimming. (Default: 20)",
    )
    parser.add_argument(
        "--max-n",
        type=float,
        default=None,
        help="Discard reads with more than this many N bases (or, if below 1.0, "
        "this proportion of Ns). Discarded reads are tagged reason=too_many_n "
        "in the discard output.",
    )
    parser.add_argument(
        "--min-avg-quality",
        type=float,
        default=None,
        help="Discard reads whose mean Phred quality is below this threshold. "
        "Discarded reads are tagged reason=low_quality in the discard output.",
    )

    parser.add_argument(
        "--with-rname-suffix",
        action="store_true",
        help="Indicate if read names have MGI-style suffixes like '/1', '/2', '.1', or '.2' to be stripped.",
    )

    parser.add_argument(
        "--ensure-inline-barcode",
        action="store_true",
        help="If set, reads without the specified inline barcode(s) will be "
        "written to the discard file, tagged reason=no_barcode. Requires "
        "adapter scheme to have inline barcodes.",
    )

    parser.add_argument(
        "--conditional-cutter",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable/disable conditional cutting for UMIs/masks. If true (default), trim UMI/mask if adapter is found or read is long. If false, UMI/mask trimming is unconditional.",
    )

    parser.add_argument(
        "--force-trim-min-length",
        type=int,
        default=50,
        help="Minimum read length to enforce UMI/mask trimming even if no adapter is found (when --conditional-cutter is true). (Default: 50)",
    )

    parser.add_argument(
        "--force-anywhere",
        action="store_true",
        help="Force adapter trimming to match anywhere in the read, not just at the ends.",
    )

    parser.add_argument(
        "-R",
        "--rename",
        "--name-format",
        dest="name_format",
        type=str,
        default=None,
        help="Custom read name template (short flag: -R). Uses cutadapt's brace "
        "variables "
        "({id}, {header}, {comment}, {cut_prefix}, {cut_suffix}, "
        "{adapter_name}, {match_sequence}, {rc}) plus positional captures "
        "{1}, {2}, ... (1-based, scheme written order; N-captures and inline "
        "barcodes counted together, e.g. YAML scheme label: names can also be "
        "used). Transform functions "
        "may wrap a capture and nest, e.g. rc({1}), upper(RC({1})), "
        "rev({1}), comp({1}), canon({1}), len({1}), left({2},6), "
        "right({2},6), slice({1},1,4). In paired mode {r1.1}/{r2.1} force a "
        "side. Example: --rename '{id}_BC1:{1}_BC2:{2}_umi:rc({3})'. "
        "--name-format is accepted as an alias.",
    )

    parser.add_argument(
        "--no-auto-inline",
        action="store_true",
        help="Disable auto-detection of inline barcodes written in uppercase "
        "between known sequencing primers. By default, an uppercase run "
        "adjacent to (or between) recognized Illumina/BGI primers is treated "
        "as an inline barcode and split/trimmed accordingly.",
    )
    parser.add_argument(
        "--list-primers",
        action="store_true",
        help="List the known sequencing primers used for inline-barcode "
        "auto-detection, then exit.",
    )

    parser.add_argument(
        "--auto-rc",
        action="store_true",
        help="Automatically reverse complement reads if the library strand "
        "(from adapter scheme) is '-'. Only applies to single-end data: for "
        "paired-end the R1/R2 orientation is fixed by the scheme's strand "
        "marker and --auto-rc is ignored (a warning is printed).",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for trimming. (Default: 1)",
    )

    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Print the sequence of modifier steps instead of running the pipeline. Does not create output files.",
    )

    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "--list-adapters",
        action="store_true",
        help="List all built-in adapter names and their schemes, then exit.",
    )

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(0)

    # Work around argparse's nargs='+' greediness (see _reorder_inputs_first):
    # allow input files to be given before *or* after -o/-d.
    args = parser.parse_args(_reorder_inputs_first(sys.argv[1:],
                                                   _option_arities(parser)))

    if args.list_adapters:
        print_builtin_adapters()
        sys.exit(0)

    if args.list_primers:
        _print_known_primers()
        sys.exit(0)

    # Check if input file is provided
    if not args.input_file:
        logging.error("Input file is required.")
        sys.exit(1)
    elif len(args.input_file) > 2:
        logging.error("Input file can not be more than two.")
        sys.exit(1)

    # Adapter resolution: -A is name-first, -a is a raw scheme.
    # A single -A/-a value is either a built-in adapter name or a custom
    # grammar scheme: a value that matches a built-in resolves to that
    # scheme; any other value is treated as the grammar scheme itself.
    if args.adapter_scheme is None:
        logging.error("Adapter scheme is required. Use -A.")
        sys.exit(1)
    builtin = BUILDIN_ADAPTERS.get(args.adapter_scheme.upper())
    if builtin is not None:
        logging.info(
            f"Resolved built-in adapter name '{args.adapter_scheme}' to scheme."
        )
        args.adapter_scheme = builtin
    else:
        logging.info(
            f"'{args.adapter_scheme}' is not a built-in adapter name; "
            f"interpreting as grammar scheme."
        )

    def validate_output_file(output_files, input_files, output_prefix, output_suffix):
        """Helper function to determine output file names."""
        default_format = ".fastq.gz"
        r1_suffix = "_" + output_suffix + "_R1" + default_format
        r2_suffix = "_" + output_suffix + "_R2" + default_format

        if output_files:  # User provided output files explicitly
            if len(output_files) != len(input_files):
                logging.error(
                    f"Number of {output_suffix} output files ({len(output_files)}) must match number of input files ({len(input_files)})."
                )
                sys.exit(1)
            return output_files
        elif output_prefix is not None:  # User provided a prefix
            if len(input_files) == 1:
                return [output_prefix + r1_suffix]
            else:
                return [output_prefix + r1_suffix, output_prefix + r2_suffix]
        else:  # Derive from input file names
            if len(input_files) == 1:
                return [remove_fq_suffix(input_files[0]) + r1_suffix]
            else:
                return [
                    remove_fq_suffix(input_files[0]) + r1_suffix,
                    remove_fq_suffix(input_files[1]) + r2_suffix,
                ]
        # This line should not be reached given the logic above.
        # However, to satisfy linters or for extreme edge cases:
        return [None] * len(input_files)

    args.output_file = validate_output_file(
        args.output_file, args.input_file, args.output_prefix, "trimmed"
    )
    args.discard_file = validate_output_file(
        args.discard_file, args.input_file, args.output_prefix, "discard"
    )
    if args.ensure_inline_barcode and args.discard_file[0] is None:
        logging.warning(
            "--ensure-inline-barcode is set but no discard output is configured "
            "(use -d/--discard-file or -O/--output-prefix). Reads that miss an "
            "inline barcode will NOT be screened out."
        )

    try:
        run_cutseq(args)
    except (ValueError, OSError) as e:
        # Scheme / --rename template / YAML / file-not-found errors: report
        # cleanly, no traceback.
        logging.error(str(e))
        sys.exit(2)


if __name__ == "__main__":
    main()
