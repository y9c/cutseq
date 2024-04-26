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
import re
import sys

from cutadapt.adapters import (
    BackAdapter,
    NonInternalBackAdapter,
    NonInternalFrontAdapter,
    PrefixAdapter,
    RightmostFrontAdapter,
    SuffixAdapter,
)
from cutadapt.files import InputPaths, OutputFiles
from cutadapt.modifiers import (
    AdapterCutter,
    PairedEndRenamer,
    QualityTrimmer,
    Renamer,
    SingleEndModifier,
    SuffixRemover,
    UnconditionalCutter,
)
from cutadapt.pipeline import PairedEndPipeline, SingleEndPipeline
from cutadapt.predicates import Predicate, TooShort
from cutadapt.report import Statistics, minimal_report
from cutadapt.runners import make_runner
from cutadapt.steps import (
    PairedEndFilter,
    PairedEndSink,
    SingleEndFilter,
    SingleEndSink,
)
from cutadapt.utils import Progress

#  monkey patching ....
original_method = Statistics._collect_modifier


def patched_problematic_method(self, *args, **kwargs):
    try:
        return original_method(self, *args, **kwargs)
    except AssertionError:
        pass


Statistics._collect_modifier = patched_problematic_method

## Enhnace function for cutadapt


class IsUntrimmedAny(Predicate):
    """
    Select reads for which no adapter match was found
    """

    def __init__(self, ref_adapters):
        self.ref_adapters = ref_adapters

    def __repr__(self):
        return "IsUntrimmedAny()"

    def test(self, read, info):
        # check not all adapters with ref_adapters are exist in the matches
        match_adapters = [match.adapter for match in info.matches]
        if any([adapter not in match_adapters for adapter in self.ref_adapters]):
            return True
        return False


class ReverseComplementConverter(SingleEndModifier):
    def __init__(self):
        self.rcd = False

    def __repr__(self):
        return "ReverseComplementConverter()"

    def __call__(self, read, info):
        return read.reverse_complement()


__version__ = importlib.metadata.version(__package__ or __name__)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s -  %(levelname)s - %(message)s",
)


def reverse_complement(b):
    return "".join(
        [dict(zip("ATGCNatgcn", "TACGNtacgn"))[x] for x in b[::-1] if x in "ATGCNatgcn"]
    )


def remove_fq_suffix(f):
    suffixes_base = ["_R1_001", "_R2_001", "_R1", "_R2", ""]
    suffixes = [
        y + "." + x for x in ["fastq.gz", "fq.gz", "fastq", "fq"] for y in suffixes_base
    ]
    print(suffixes)
    for suffix in suffixes:
        if f.endswith(suffix):
            return f.removesuffix(suffix)
    return f


class BarcodeSeq:
    def __init__(self, seq):
        self.fw = seq
        self.rc = reverse_complement(seq)
        self.len = len(seq)


class BarcodeConfig:
    """
    Adapter scheme:
    (p5)(inline5)(umi5)(mask5)(strand)(mask3)(umi3)(inline3)(p7)
    """

    def __init__(self, adapter=None):
        self.strand = None
        self.p5 = BarcodeSeq("")
        self.p7 = BarcodeSeq("")
        self.inline5 = BarcodeSeq("")
        self.inline3 = BarcodeSeq("")
        self.umi5 = BarcodeSeq("")
        self.umi3 = BarcodeSeq("")
        self.mask5 = BarcodeSeq("")
        self.mask3 = BarcodeSeq("")
        if adapter is not None:
            self._parse_barcode(adapter)

    def _parse_barcode(self, b):
        m = re.match(
            r"(?P<p5>[ATGCatgc]+)(\((?P<inline5>[ATGCatgc]+)\))?(?P<umi5>N*)(?P<mask5>X*)(?P<strand>-|>|<)(?P<mask3>X*)(?P<umi3>N*)(\((?P<inline3>[ATGCatgc]+)\))?(?P<p7>[ATGCatgc]+)",
            b,
        )
        if m is None:
            logging.error(f"barcode {b} is not valid")
            sys.exit(1)
        d = m.groupdict()
        if d["inline5"] is None:
            d["inline5"] = ""
        if d["inline3"] is None:
            d["inline3"] = ""
        self.strand = "+" if d["strand"] == ">" else "-" if d["strand"] == "<" else None
        self.p5 = BarcodeSeq(d["p5"])
        self.p7 = BarcodeSeq(d["p7"])
        self.inline5 = BarcodeSeq(d["inline5"])
        self.inline3 = BarcodeSeq(d["inline3"])
        self.umi5 = BarcodeSeq(d["umi5"])
        self.umi3 = BarcodeSeq(d["umi3"])
        self.mask5 = BarcodeSeq(d["mask5"])
        self.mask3 = BarcodeSeq(d["mask3"])


class CutadaptConfig:
    def __init__(self):
        self.rname_suffix = False
        self.ensure_inline_barcode = False
        self.trim_polyA = False
        self.min_length = 20
        self.min_quality = 20
        self.auto_rc = False
        self.dry_run = False
        self.threads = 1
        self.json_file = None


BUILDIN_ADAPTERS = {
    "INLINE": "AGTTCTACAGTCCGACGATCNNNNN>NNNNN(ATCACG)AGATCGGAAGAGCACACGTC",
    # Small RNA, double ligation method, without barcode
    # p5 - insert - p7
    # (Optional) trim 2nt on both end to increase quality
    "SMALLRNA": "CACGACGCTCTTCCGATCT>AGATCGGAAGAGCACACGTC",
    # p5 - (random rt tail in TSO) - reverse insert - (random primer start?) - p7
    "TAKARAV2": "ACACGACGCTCTTCCGATCTXXX<XXXAGATCGGAAGAGCACACGTC",
    # p5 - (random rt tail in ligation) - reverse insert - (random primer start?) - p7
    "STRANDED": "ACACGACGCTCTTCCGATCTX<XXXAGATCGGAAGAGCACACGTC",
    # p5 - reverse insert - 14ntUMI - p7
    # 14nt UMI = (8 nt UMIs + 3 nt UMI linker + 3 nt from Pico v3 SMART UMI Adapter)
    # IMPORTANT: The UMI liker and UMI adapter can be different, even the 8nt UMI is the same. very weired.
    # NOTE: if insert is too short, also need to add -u -14 to trim readthrough in R1
    "TAKARAV3": "ACACGACGCTCTTCCGATCTXXX<XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC",
    # p5 - [might be 6bp of polyC] - reverse insert (cDNA) - adaptase tail (CCCCCC) - p7
    # 6nt of polyG in 5' of R1 might from random RT primer
    # adaptase tail can be as long as 15bp at the 5' of R2 of polyG)
    # no UMI, but try to use random polyC tail as UMI
    "SWIFT": "ACACGACGCTCTTCCGATCTXXXXXX<XXXXXXXXXXXXXXXAGATCGGAAGAGCACACGTC",
}


def pipeline_single(input1, output1, short1, untrimmed1, barcode, settings):
    modifiers = []
    # step 1: remove suffix in the read name
    modifiers.extend([SuffixRemover(".1"), SuffixRemover("/1")])
    # step 2: remove adapter on the 5' end, artifact of template switching
    modifiers.append(
        AdapterCutter(
            [
                RightmostFrontAdapter(
                    sequence=barcode.p5.fw, max_errors=0.25, min_overlap=10
                )
            ],
            times=1,
        )
    )
    # step 3: remove adapter on the 3' end, read though in the sequencing
    modifiers.append(
        AdapterCutter(
            [BackAdapter(sequence=barcode.p7.fw, max_errors=0.2, min_overlap=3)],
            times=2,
        ),
    )
    # step 4: trim inline barcode
    if barcode.inline5.len > 0:
        adapter_inline5 = PrefixAdapter(sequence=barcode.inline5.fw, max_errors=0.2)
        modifiers.append(AdapterCutter([adapter_inline5], times=1))
    else:
        adapter_inline5 = None
    if barcode.inline3.len > 0:
        adapter_inline3 = SuffixAdapter(sequence=barcode.inline3.fw, max_errors=0.2)
        modifiers.append(AdapterCutter([adapter_inline3], times=1))
    else:
        adapter_inline3 = None

    # step 5: extract UMI
    if barcode.umi5.len > 0:
        modifiers.append(UnconditionalCutter(barcode.umi5.len))
    if barcode.umi3.len > 0:
        modifiers.append(UnconditionalCutter(-barcode.umi3.len))
    if barcode.umi5.len + barcode.umi3.len > 0:
        modifiers.append(Renamer("{id}_{cut_prefix}{cut_suffix}"))
    else:
        modifiers.append(Renamer("{id}"))

    # step 6: mask tail in the RNA, which might be artifact of RT
    if barcode.mask5.len > 0:
        modifiers.append(UnconditionalCutter(barcode.mask5.len))
    if barcode.mask3.len > 0:
        modifiers.append(UnconditionalCutter(-barcode.mask3.len))
    # step 7: trim polyA
    if settings.trim_polyA:
        if barcode.strand == "+":
            modifiers.append(
                AdapterCutter(
                    [NonInternalBackAdapter(sequence="A" * 100, max_errors=0.15)]
                )
            )
        elif barcode.strand == "-":
            modifiers.append(
                AdapterCutter(
                    [NonInternalFrontAdapter(sequence="T" * 100, max_errors=0.15)]
                )
            )
        else:
            logging.info("No strand information provided, skip polyA trimming.")
    # step 8: quality control, remove short reads
    modifiers.append(
        QualityTrimmer(cutoff_front=0, cutoff_back=settings.min_quality),
    )

    # step 9: reverse complement the read
    if settings.auto_rc:
        if barcode.strand == "-":
            modifiers.append(ReverseComplementConverter())
        else:
            logging.warn(
                "Library is not (-) strand, but --auto-rc is enabled. Ignored."
            )

    # dry run and exit code
    if settings.dry_run:
        for i, m in enumerate(modifiers, 1):
            print(f"Step {i}: {m}")
        return

    inpaths = InputPaths(input1)

    with make_runner(inpaths, cores=settings.threads) as runner:
        outfiles = OutputFiles(
            proxied=settings.threads > 1,
            qualities=runner.input_file_format().has_qualities(),
            interleaved=False,
        )
        steps = []
        # TODO: report info
        # --info-file=info.txt
        # PairedSingleEndStep(InfoFileWriter(outfiles.open_text("info.txt"))),
        steps.append(
            # -m 10
            SingleEndFilter(
                TooShort(settings.min_length), outfiles.open_record_writer(short1)
            ),
        )
        # TODO: --max-n=0 support
        if (
            barcode.inline5.len + barcode.inline3.len > 0
            and settings.ensure_inline_barcode
        ) or (untrimmed1 is not None):
            ref_adapters = []
            if adapter_inline5 is not None:
                ref_adapters.append(adapter_inline5)
            if adapter_inline3 is not None:
                ref_adapters.append(adapter_inline3)
            steps.append(
                SingleEndFilter(
                    IsUntrimmedAny(ref_adapters),
                    outfiles.open_record_writer(untrimmed1, interleaved=False),
                ),
            )
        steps.append(
            # -o
            SingleEndSink(outfiles.open_record_writer(output1, interleaved=False)),
        )
        pipeline = SingleEndPipeline(modifiers, steps)
        stats = runner.run(pipeline, Progress(), outfiles)
        if settings.json_file:
            with open(settings.json_file, "w") as json_file:
                json_file.write(json.dumps(stats.as_json(), indent=2))
        print(minimal_report(stats, time=None, gc_content=None), file=sys.stderr)
    outfiles.close()


def pipeline_paired(
    input1,
    input2,
    output1,
    output2,
    short1,
    short2,
    untrimmed1,
    untrimmed2,
    barcode,
    settings,
):
    modifiers = []
    # step 1: remove suffix in the read name
    modifiers.extend(
        [
            (SuffixRemover(".1"), SuffixRemover(".2")),
            (SuffixRemover("/1"), SuffixRemover("/2")),
        ]
    )
    # step 2: remove adapter on the 5' end, artifact of template switching
    modifiers.append(
        (
            AdapterCutter(
                [
                    RightmostFrontAdapter(
                        sequence=barcode.p5.fw, max_errors=0.25, min_overlap=10
                    )
                ],
                times=1,
            ),
            AdapterCutter(
                [
                    RightmostFrontAdapter(
                        sequence=barcode.p7.rc, max_errors=0.25, min_overlap=10
                    )
                ],
                times=1,
            ),
        ),
    )
    # step 3: remove adapter on the 3' end, read though in the sequencing
    modifiers.append(
        (
            AdapterCutter(
                [BackAdapter(sequence=barcode.p7.fw, max_errors=0.2, min_overlap=3)],
                times=2,
            ),
            AdapterCutter(
                [BackAdapter(sequence=barcode.p5.rc, max_errors=0.2, min_overlap=3)],
                times=2,
            ),
        ),
    )
    # step 4: trim inline barcode
    if barcode.inline5.len > 0:
        adapter_inline5 = PrefixAdapter(sequence=barcode.inline5.fw, max_errors=0.2)
        modifiers.append(
            (
                AdapterCutter([adapter_inline5], times=1),
                UnconditionalCutter(-barcode.inline5.len),
            )
        )
    else:
        adapter_inline5 = None
    if barcode.inline3.len > 0:
        adapter_inline3 = PrefixAdapter(sequence=barcode.inline3.rc, max_errors=0.2)
        modifiers.append(
            (
                UnconditionalCutter(-barcode.inline3.len),
                AdapterCutter([adapter_inline3], times=1),
            )
        )
    else:
        adapter_inline3 = None

    # step 5: extract UMI
    if barcode.umi5.len > 0:
        modifiers.append(
            (
                UnconditionalCutter(barcode.umi5.len),
                UnconditionalCutter(-barcode.umi5.len),
            ),
        )
    if barcode.umi3.len > 0:
        modifiers.append(
            (
                UnconditionalCutter(-barcode.umi3.len),
                UnconditionalCutter(barcode.umi3.len),
            )
        )
    if barcode.umi5.len + barcode.umi3.len > 0:
        modifiers.append(PairedEndRenamer("{id}_{r1.cut_prefix}{r2.cut_prefix}"))
    else:
        modifiers.append(PairedEndRenamer("{id}"))

    # step 6: mask tail in the RNA, which might be artifact of RT
    if barcode.mask5.len > 0:
        modifiers.append(
            (
                UnconditionalCutter(barcode.mask5.len),
                UnconditionalCutter(-barcode.mask5.len),
            )
        )
    if barcode.mask3.len > 0:
        modifiers.append(
            (
                UnconditionalCutter(-barcode.mask3.len),
                UnconditionalCutter(barcode.mask3.len),
            )
        )
    # step 7: trim polyA
    if settings.trim_polyA:
        if barcode.strand == "+":
            modifiers.append(
                (
                    AdapterCutter(
                        [NonInternalBackAdapter(sequence="A" * 100, max_errors=0.15)]
                    ),
                    AdapterCutter(
                        [NonInternalFrontAdapter(sequence="T" * 100, max_errors=0.15)]
                    ),
                )
            )
        elif barcode.strand == "-":
            modifiers.append(
                (
                    AdapterCutter(
                        [NonInternalFrontAdapter(sequence="T" * 100, max_errors=0.15)]
                    ),
                    AdapterCutter(
                        [NonInternalBackAdapter(sequence="A" * 100, max_errors=0.15)]
                    ),
                )
            )
        else:
            logging.info("No strand information provided, skip polyA trimming.")
    # step 8: quality control, remove short reads
    modifiers.append(
        (
            QualityTrimmer(cutoff_front=0, cutoff_back=settings.min_quality),
            QualityTrimmer(cutoff_front=0, cutoff_back=settings.min_quality),
        )
    )

    # step 9: reverse complement the read
    # NOTE: Do not need to rc, switch R1 and R2 since it is paired-end data
    if settings.auto_rc:
        if barcode.strand != "-":
            logging.warn(
                "Library is not (-) strand, but --auto-rc is enabled. Ignored."
            )

    # dry run and exit code
    if settings.dry_run:
        for i, m in enumerate(modifiers, 1):
            print(f"Step {i}: {m}")
        return

    inpaths = InputPaths(input1, input2)

    with make_runner(inpaths, cores=settings.threads) as runner:
        outfiles = OutputFiles(
            proxied=settings.threads > 1,
            qualities=runner.input_file_format().has_qualities(),
            interleaved=False,
        )
        steps = []
        # --info-file=info.txt
        # PairedSingleEndStep(InfoFileWriter(outfiles.open_text("info.txt"))),
        # -m 10:10
        steps.append(
            PairedEndFilter(
                TooShort(settings.min_length),
                TooShort(settings.min_length),
                outfiles.open_record_writer(short1, short2, interleaved=False),
            )
        )
        # TODO: --max-n=0 support
        if (
            barcode.inline5.len + barcode.inline3.len > 0
            and settings.ensure_inline_barcode
        ) or (untrimmed1 is not None and untrimmed2 is not None):
            steps.append(
                PairedEndFilter(
                    IsUntrimmedAny([adapter_inline5] if adapter_inline5 else []),
                    IsUntrimmedAny([adapter_inline3] if adapter_inline3 else []),
                    outfiles.open_record_writer(
                        untrimmed1, untrimmed2, interleaved=False
                    ),
                    pair_filter_mode="any",
                )
            )
        steps.append(
            # -o ... -p ...
            PairedEndSink(
                outfiles.open_record_writer(output1, output2)
                if (settings.auto_rc and barcode.strand == "-")
                else outfiles.open_record_writer(output2, output1)
            )
        )
        pipeline = PairedEndPipeline(modifiers, steps)
        stats = runner.run(pipeline, Progress(), outfiles)
        if settings.json_file:
            with open(settings.json_file, "w") as json_file:
                json_file.write(json.dumps(stats.as_json(), indent=2))
        print(minimal_report(stats, time=None, gc_content=None), file=sys.stderr)

    outfiles.close()


def run_cutseq(args):
    barcode_config = BarcodeConfig(args.adapter_scheme)
    settings = CutadaptConfig()
    settings.rname_suffix = args.with_rname_suffix
    settings.ensure_inline_barcode = args.ensure_inline_barcode
    settings.trim_polyA = args.trim_polyA
    settings.threads = args.threads
    settings.min_length = args.min_length
    settings.dry_run = args.dry_run
    settings.auto_rc = args.auto_rc
    settings.json_file = args.json_file
    if len(args.input_file) == 1:
        pipeline_single(
            args.input_file[0],
            args.output_file[0],
            args.short_file[0],
            args.untrimmed_file[0],
            barcode_config,
            settings,
        )
    else:
        pipeline_paired(
            args.input_file[0],
            args.input_file[1],
            args.output_file[0],
            args.output_file[1],
            args.short_file[0],
            args.short_file[1],
            args.untrimmed_file[0],
            args.untrimmed_file[1],
            barcode_config,
            settings,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Trim sequencing adapters from NGS data automatically."
    )
    # input file can be one or two for single or paired-end reads, but can not be more than two
    parser.add_argument(
        "input_file",
        type=str,
        nargs="+",
        help="Input file path for NGS data, one or two files.",
    )
    # output file can be number of files matching the input files, if not provided it will generate based on the output prefix,
    # if no output prefix provided it will generate based on the input file name
    parser.add_argument(
        "-a",
        "--adapter-scheme",
        type=str,
        help="Adapter sequence configuration.",
    )
    parser.add_argument(
        "-A",
        "--adapter-name",
        type=str.upper,
        choices=BUILDIN_ADAPTERS.keys(),
        help="Built-in adapter name.",
    )
    parser.add_argument(
        "-O",
        "--output-prefix",
        type=str,
        help="Output file prefix for keep trimmed data.",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        nargs="+",
        help="Output file path for keep trimmed data.",
    )
    # discard short reads
    parser.add_argument(
        "-s",
        "--short-file",
        type=str,
        nargs="+",
        help="Output file path for discarded too short data.",
    )
    # discard inline barcode untrimmed reads
    parser.add_argument(
        "-u",
        "--untrimmed-file",
        type=str,
        nargs="+",
        help="Output file path for discarded reads without inline barcode.",
    )
    parser.add_argument(
        "--json-file",
        type=str,
        help="Output json file for statistics.",
    )

    parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=20,
        help="Minimum quality of the read tails in the reads to keep.",
    )
    parser.add_argument(
        "-m",
        "--min-length",
        type=int,
        default=20,
        help="Minimum length of the reads to keep.",
    )

    parser.add_argument(
        "--with-rname-suffix",
        action="store_true",
        help="R1 and R2 suffix cotains suffix. MGI platform.",
    )

    parser.add_argument(
        "--ensure-inline-barcode",
        action="store_true",
        help="Output discarded reads without inline barcode.",
    )

    parser.add_argument("--trim-polyA", action="store_true", help="Trim polyA tail.")

    parser.add_argument(
        "--auto-rc",
        action="store_true",
        help="Reverse complement (-) strand automatically if the library is reverse orientation.",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for trimming.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Print command instead of running it.",
    )

    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(0)

    args = parser.parse_args()

    if args.adapter_name is not None:
        if args.adapter_scheme is not None:
            logging.info("Adapter scheme is provided, ignore adapter name.")
        else:
            args.adapter_scheme = BUILDIN_ADAPTERS.get(args.adapter_name.upper())
            if args.adapter_scheme is None:
                logging.error("Adapter name is not valid.")
                sys.exit(1)
    elif args.adapter_scheme is None:
        logging.error("Adapter scheme or name is required.")
        sys.exit(1)
    args.adapter_scheme = args.adapter_scheme.replace(" ", "").upper()

    if len(args.input_file) > 2:
        logging.error("Input file can not be more than two.")
        sys.exit(1)

    def validate_output_file(output_files, input_files, output_prefix, output_suffix):
        default_format = ".fastq.gz"
        r1 = "_" + output_suffix + "_R1" + default_format
        r2 = "_" + output_suffix + "_R2" + default_format
        if output_files:
            if len(output_files) != len(input_files):
                logging.error("Output file should be same as input file.")
                sys.exit(1)
            return output_files
        elif output_prefix is not None:
            if len(input_files) == 1:
                return [output_prefix + r1]
            else:
                return [output_prefix + r1, output_prefix + r2]
        else:
            if len(input_files) == 1:
                return [remove_fq_suffix(input_files[0]) + r1]
            else:
                return [
                    remove_fq_suffix(input_files[0]) + r1,
                    remove_fq_suffix(input_files[1]) + r2,
                ]
        return output_files

    args.output_file = validate_output_file(
        args.output_file, args.input_file, args.output_prefix, "trimmed"
    )
    args.short_file = validate_output_file(
        args.short_file, args.input_file, args.output_prefix, "short"
    )

    def _check_with_inline_barcode(s):
        # inline barcode is in bracket () and length > 0
        return re.match(r".*\([ATGCatgc]+\).*", s) is not None

    args.untrimmed_file = (
        validate_output_file(
            args.untrimmed_file, args.input_file, args.output_prefix, "untrimmed"
        )
        if (
            args.untrimmed_file is not None
            or (
                _check_with_inline_barcode(args.adapter_scheme)
                and args.ensure_inline_barcode
            )
        )
        else [None, None]
    )

    run_cutseq(args)


if __name__ == "__main__":
    main()
