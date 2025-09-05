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

import cutadapt
from cutadapt.adapters import (
    BackAdapter,
    NonInternalBackAdapter,
    NonInternalFrontAdapter,
    PrefixAdapter,
    RightmostFrontAdapter,
    SuffixAdapter,
)
from cutadapt.files import InputPaths, OutputFiles
from cutadapt.info import ModificationInfo
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

from .common import (
    BarcodeConfig,
    load_adapters,
    print_builtin_adapters,
    remove_fq_suffix,
)

# Initialize built-in adapters at the top so it's available for argument parsing
BUILDIN_ADAPTERS = load_adapters()

#  monkey patching ....
original_method = Statistics._collect_modifier


def patched_problematic_method(self, *args, **kwargs):
    """
    Patches Statistics._collect_modifier to ignore AssertionError.
    This is a workaround for a potential issue in cutadapt's statistics collection.
    """
    try:
        return original_method(self, *args, **kwargs)
    except AssertionError:
        pass


Statistics._collect_modifier = patched_problematic_method

## Enhnace function for cutadapt


class IsUntrimmedAny(Predicate):
    """
    Cutadapt Predicate to select reads where at least one of the specified
    reference adapters was *not* found among the read's matches.

    This is useful for filtering reads that are expected to have certain inline
    barcodes, allowing reads that missed one or more of these barcodes to be
    written to a separate file.

    :param ref_adapters: A list of cutadapt adapter objects that are expected to be found.
    :type ref_adapters: list[cutadapt.adapters.Adapter]
    """

    def __init__(self, ref_adapters):
        self.ref_adapters = ref_adapters

    def __repr__(self):
        return f"IsUntrimmedAny(ref_adapters={self.ref_adapters!r})"

    def test(self, read, info):
        """
        Tests if any of the reference adapters are missing from the read's matches.

        :param read: The read object.
        :param info: The ModificationInfo object for the read.
        :return: True if at least one reference adapter is not in info.matches, False otherwise.
        :rtype: bool
        """
        # check not all adapters with ref_adapters are exist in the matches
        match_adapters = [match.adapter for match in info.matches]
        if any([adapter not in match_adapters for adapter in self.ref_adapters]):
            return True
        return False


class ConditionalCutter(SingleEndModifier):
    """
    Cutadapt SingleEndModifier that conditionally cuts a fixed length from
    the start or end of a read.

    The primary condition is based on whether adapters were matched (`info.matches`).
    If no adapters were matched AND the read is shorter than `force_trim_min_length`,
    the read is returned unmodified.
    Otherwise (if adapters were matched OR the read is long enough), the specified
    `length` is cut.
    - If `length` is positive, it's cut from the 5' end.
    - If `length` is negative, it's cut from the 3' end.

    This modifier is typically used for UMI or mask trimming where trimming should
    occur if an adapter was found, or if the read is long enough to suggest that
    the UMI/mask is present even without an adapter match (e.g., due to sequencing
    errors in the adapter region).

    :param length: The length of the sequence to cut. Positive for 5' end, negative for 3' end.
    :type length: int
    :param force_trim_min_length: The minimum read length to enforce trimming even
                                  if no adapter was found.
    :type force_trim_min_length: int
    """

    def __init__(self, length: int, force_trim_min_length: int = 50):
        self.length = length
        self.force_trim_min_length = force_trim_min_length

    def __repr__(self):
        return f"ConditionalCutter(length={self.length}, force_trim_min_length={self.force_trim_min_length})"

    def __call__(self, read, info: ModificationInfo):
        """
        Applies the conditional cut to the read.

        :param read: The read object.
        :param info: The ModificationInfo object for the read.
        :return: The modified read.
        :rtype: cutadapt.sequence.Sequence
        """
        if not info.matches and len(read.sequence) < self.force_trim_min_length:
            return read
        if self.length > 0:
            info.cut_prefix = read.sequence[: self.length]
            return read[self.length :]
        elif self.length < 0:
            info.cut_suffix = read.sequence[self.length :]
            return read[: self.length]


class ReverseComplementConverter(SingleEndModifier):
    """
    Cutadapt SingleEndModifier that reverse complements a read.

    This modifier is used to orient reads consistently, typically when the
    library preparation protocol results in reads from the '-' strand.
    """

    def __init__(self):
        self.rcd = False  # This attribute seems unused but kept for compatibility if needed later.

    def __repr__(self):
        return "ReverseComplementConverter()"

    def __call__(self, read, info):
        """
        Reverse complements the input read.

        :param read: The read object.
        :param info: The ModificationInfo object for the read (unused by this modifier).
        :return: The reverse complemented read.
        :rtype: cutadapt.sequence.Sequence
        """
        return read.reverse_complement()


__version__ = importlib.metadata.version(__package__ or __name__)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s -  %(levelname)s - %(message)s",
)


class CutadaptConfig:
    """
    Configuration class to store settings for the cutadapt pipeline.

    This class centralizes various options that control the behavior of
    the adapter trimming and filtering process.
    """

    def __init__(self):
        self.rname_suffix = False
        self.ensure_inline_barcode = False
        self.trim_polyA = False
        self.trim_polyA_wo_direction = False
        self.conditional_cutter = True
        self.min_length = 20
        self.min_quality = 20
        self.auto_rc = False
        self.dry_run = False
        self.threads = 1
        self.json_file = None
        self.force_trim_min_length = 50  # Default value
        self.force_anywhere = False


def json_report(
    file,
    stats,
    barcode,
    input1,
    input2,
    output1,
    output2,
    short1,
    short2,
    untrimmed1,
    untrimmed2,
):
    """
    Generates a JSON report summarizing the trimming statistics.

    The report includes cutadapt version, input/output file paths, barcode configuration,
    and detailed trimming statistics from the cutadapt run.

    :param file: Path to the output JSON file.
    :type file: str
    :param stats: Cutadapt Statistics object.
    :type stats: cutadapt.report.Statistics
    :param barcode: BarcodeConfig object used for the run.
    :type barcode: cutseq.common.BarcodeConfig
    :param input1: Path to the first input FASTQ file.
    :type input1: str
    :param input2: Path to the second input FASTQ file (None for single-end).
    :type input2: str, optional
    :param output1: Path to the first output trimmed FASTQ file.
    :type output1: str
    :param output2: Path to the second output trimmed FASTQ file (None for single-end).
    :type output2: str, optional
    :param short1: Path to the first output file for short reads (None if not used).
    :type short1: str, optional
    :param short2: Path to the second output file for short reads (None if not used).
    :type short2: str, optional
    :param untrimmed1: Path to the first output file for untrimmed reads (None if not used).
    :type untrimmed1: str, optional
    :param untrimmed2: Path to the second output file for untrimmed reads (None if not used).
    :type untrimmed2: str, optional
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
            "short1": short1,
            "short2": short2,
            "untrimmed1": untrimmed1,
            "untrimmed2": untrimmed2,
        },
        "barcode": barcode.to_dict(),
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


def pipeline_single(input1, output1, short1, untrimmed1, barcode, settings):
    """
    Configures and runs the cutadapt pipeline for single-end reads.

    This function defines a series of modification and filtering steps
    based on the provided barcode configuration and settings.

    :param input1: Path to the input FASTQ file.
    :type input1: str
    :param output1: Path to the output trimmed FASTQ file.
    :type output1: str
    :param short1: Path to the output file for reads that are too short.
    :type short1: str, optional
    :param untrimmed1: Path to the output file for reads that were not trimmed
                       (e.g., missing inline barcodes).
    :type untrimmed1: str, optional
    :param barcode: BarcodeConfig object detailing the adapter scheme.
    :type barcode: cutseq.common.BarcodeConfig
    :param settings: CutadaptConfig object with pipeline settings.
    :type settings: CutadaptConfig
    """
    max_errors = 0.2
    repeat_trim_times = 1
    modifiers = []
    # step 1: remove suffix in the read name
    modifiers.extend([SuffixRemover(".1"), SuffixRemover("/1")])
    # step 2: remove adapter on the 5' end, artifact of template switching
    modifiers.append(
        AdapterCutter(
            [
                RightmostFrontAdapter(
                    sequence=barcode.p5.fw, max_errors=max_errors, min_overlap=10
                )
            ],
            times=repeat_trim_times,
        )
    )
    # step 3: remove adapter on the 3' end, read though in the sequencing
    modifiers.append(
        AdapterCutter(
            [
                BackAdapter(
                    sequence=barcode.p7.fw,
                    max_errors=max_errors,
                    min_overlap=3,
                    force_anywhere=settings.force_anywhere,
                )
            ],
            times=repeat_trim_times,
        ),
    )
    # step 4: trim inline barcode
    if barcode.inline5.len > 0:
        adapter_inline5 = PrefixAdapter(
            sequence=barcode.inline5.fw, max_errors=max_errors
        )
        modifiers.append(AdapterCutter([adapter_inline5], times=1))
    else:
        adapter_inline5 = None
    if barcode.inline3.len > 0:
        adapter_inline3 = SuffixAdapter(
            sequence=barcode.inline3.fw, max_errors=max_errors
        )
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
        pA_max_errors = 0.15
        pA_max_length = 100
        m_fwd = AdapterCutter(
            [
                NonInternalBackAdapter(
                    sequence="A" * pA_max_length, max_errors=pA_max_errors
                )
            ]
        )
        m_rev = AdapterCutter(
            [
                NonInternalFrontAdapter(
                    sequence="T" * pA_max_length, max_errors=pA_max_errors
                )
            ]
        )
        if settings.trim_polyA_wo_direction:
            modifiers.append(m_fwd)
            modifiers.append(m_rev)
        elif barcode.strand == "+":
            modifiers.append(m_fwd)
        elif barcode.strand == "-":
            modifiers.append(m_rev)
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
            logging.warning(
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

        if settings.json_file is not None:
            json_report(
                settings.json_file,
                stats,
                barcode,
                input1,
                None,
                output1,
                None,
                short1,
                None,
                untrimmed1,
                None,
            )
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
    """
    Configures and runs the cutadapt pipeline for paired-end reads.

    This function defines a series of modification and filtering steps
    for paired-end data, considering R1 and R2 reads, based on the
    provided barcode configuration and settings.

    :param input1: Path to the R1 input FASTQ file.
    :type input1: str
    :param input2: Path to the R2 input FASTQ file.
    :type input2: str
    :param output1: Path to the R1 output trimmed FASTQ file.
    :type output1: str
    :param output2: Path to the R2 output trimmed FASTQ file.
    :type output2: str
    :param short1: Path to the R1 output file for reads that are too short.
    :type short1: str, optional
    :param short2: Path to the R2 output file for reads that are too short.
    :type short2: str, optional
    :param untrimmed1: Path to the R1 output file for untrimmed reads.
    :type untrimmed1: str, optional
    :param untrimmed2: Path to the R2 output file for untrimmed reads.
    :type untrimmed2: str, optional
    :param barcode: BarcodeConfig object detailing the adapter scheme.
    :type barcode: cutseq.common.BarcodeConfig
    :param settings: CutadaptConfig object with pipeline settings.
    :type settings: CutadaptConfig
    """
    max_errors = 0.2
    repeat_trim_times = 1
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
                        sequence=barcode.p5.fw, max_errors=max_errors, min_overlap=10
                    )
                ],
                times=1,
            ),
            AdapterCutter(
                [
                    RightmostFrontAdapter(
                        sequence=barcode.p7.rc, max_errors=max_errors, min_overlap=10
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
                [
                    BackAdapter(
                        sequence=barcode.p7.fw,
                        max_errors=max_errors,
                        min_overlap=3,
                        force_anywhere=settings.force_anywhere,
                    )
                ],
                times=repeat_trim_times,
            ),
            AdapterCutter(
                [
                    BackAdapter(
                        sequence=barcode.p5.rc,
                        max_errors=max_errors,
                        min_overlap=3,
                        force_anywhere=settings.force_anywhere,
                    )
                ],
                times=repeat_trim_times,
            ),
        ),
    )
    # step 4: trim inline barcode
    if barcode.inline5.len > 0:
        adapter_inline5 = PrefixAdapter(
            sequence=barcode.inline5.fw, max_errors=max_errors
        )
        modifiers.append(
            (
                AdapterCutter([adapter_inline5], times=1),
                UnconditionalCutter(-barcode.inline5.len),
            )
        )
    else:
        adapter_inline5 = None
    if barcode.inline3.len > 0:
        adapter_inline3 = PrefixAdapter(
            sequence=barcode.inline3.rc, max_errors=max_errors
        )
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
                ConditionalCutter(
                    -barcode.umi5.len,
                    force_trim_min_length=settings.force_trim_min_length,
                )
                if settings.conditional_cutter
                else UnconditionalCutter(-barcode.umi5.len),
            ),
        )
    if barcode.umi3.len > 0:
        modifiers.append(
            (
                ConditionalCutter(
                    -barcode.umi3.len,
                    force_trim_min_length=settings.force_trim_min_length,
                )
                if settings.conditional_cutter
                else UnconditionalCutter(-barcode.umi3.len),
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
                ConditionalCutter(
                    -barcode.mask5.len,
                    force_trim_min_length=settings.force_trim_min_length,
                )
                if settings.conditional_cutter
                else UnconditionalCutter(-barcode.mask5.len),
            )
        )
    if barcode.mask3.len > 0:
        modifiers.append(
            (
                ConditionalCutter(
                    -barcode.mask3.len,
                    force_trim_min_length=settings.force_trim_min_length,
                )
                if settings.conditional_cutter
                else UnconditionalCutter(-barcode.mask3.len),
                UnconditionalCutter(barcode.mask3.len),
            )
        )
    # step 7: trim polyA
    if settings.trim_polyA:
        pA_max_errors = 0.15
        pA_max_length = 100
        m_fwd = (
            AdapterCutter(
                [
                    NonInternalBackAdapter(
                        sequence="A" * pA_max_length, max_errors=pA_max_errors
                    )
                ]
            ),
            AdapterCutter(
                [
                    NonInternalFrontAdapter(
                        sequence="T" * pA_max_length, max_errors=pA_max_errors
                    )
                ]
            ),
        )
        m_rev = (
            AdapterCutter(
                [
                    NonInternalFrontAdapter(
                        sequence="T" * pA_max_length, max_errors=pA_max_errors
                    )
                ]
            ),
            AdapterCutter(
                [
                    NonInternalBackAdapter(
                        sequence="A" * pA_max_length, max_errors=pA_max_errors
                    )
                ]
            ),
        )
        if settings.trim_polyA_wo_direction:
            modifiers.append(m_fwd)
            modifiers.append(m_rev)
        elif barcode.strand == "+":
            modifiers.append(m_fwd)
        elif barcode.strand == "-":
            modifiers.append(m_rev)
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
            logging.warning(
                "Library is not (-) strand, but --auto-rc is enabled. Ignored."
            )

    # dry run and exit code
    if settings.dry_run:
        for b in [
            "p5",
            "p7",
            "inline5",
            "inline3",
            "umi5",
            "umi3",
            "mask5",
            "mask3",
            "strand",
        ]:
            print(f"{b}: {getattr(barcode, b)}")
        for i, m in enumerate(modifiers, 1):
            logging.info(f"Step {i}: {m}")
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
                outfiles.open_record_writer(output2, output1)
                if (settings.auto_rc and barcode.strand == "-")
                else outfiles.open_record_writer(output1, output2)
            )
        )
        pipeline = PairedEndPipeline(modifiers, steps)
        stats = runner.run(pipeline, Progress(), outfiles)

        if settings.json_file is not None:
            json_report(
                settings.json_file,
                stats,
                barcode,
                input1,
                input2,
                output1,
                output2,
                short1,
                short2,
                untrimmed1,
                untrimmed2,
            )
        print(minimal_report(stats, time=None, gc_content=None), file=sys.stderr)

    outfiles.close()


def run_cutseq(args):
    """
    Sets up configurations and executes the appropriate cutadapt pipeline.

    This function initializes the BarcodeConfig and CutadaptConfig based on
    the command-line arguments. It then determines whether to run the
    single-end or paired-end pipeline based on the number of input files
    provided.

    :param args: Parsed command-line arguments from argparse.
    :type args: argparse.Namespace
    """
    barcode_config = BarcodeConfig(args.adapter_scheme)
    settings = CutadaptConfig()
    settings.rname_suffix = args.with_rname_suffix
    settings.ensure_inline_barcode = args.ensure_inline_barcode
    settings.trim_polyA = args.trim_polyA
    settings.trim_polyA_wo_direction = args.trim_polyA_wo_direction
    settings.conditional_cutter = args.conditional_cutter
    settings.threads = args.threads
    settings.min_length = args.min_length
    settings.min_quality = args.min_quality
    settings.dry_run = args.dry_run
    settings.auto_rc = args.auto_rc
    settings.json_file = args.json_file
    settings.force_trim_min_length = args.force_trim_min_length  # Pass from args
    settings.force_anywhere = args.force_anywhere
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
    """
    Parses command-line arguments and initiates the adapter trimming process.

    This is the main entry point for the `cutseq run` script. It sets up
    the argument parser, validates arguments, prepares output file names,
    and then calls `run_cutseq` to perform the trimming.
    """
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
        "-a",
        "--adapter-scheme",
        type=str,
        help="Adapter sequence configuration string. Example: P5(INLINE5)UMI5XXXS>P7(INLINE3)UMI3XXXS. "
        "Where P5/P7 are adapter sequences, (INLINE5/3) are optional inline barcodes, "
        "UMI5/3 are N's for UMI bases, XXX are mask sequences, S is strand (>/< or -).",
    )
    parser.add_argument(
        "-A",
        "--adapter-name",
        help="Built-in adapter name. choices:\n" + ",".join(BUILDIN_ADAPTERS.keys()),
    )
    parser.add_argument(
        "-O",
        "--output-prefix",
        type=str,
        help="Output file prefix for trimmed, short, and untrimmed data. "
        "If not provided, output filenames are derived from input filenames.",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        nargs="+",
        help="Output file path(s) for successfully trimmed reads. Must match number of input files.",
    )
    # discard short reads
    parser.add_argument(
        "-s",
        "--short-file",
        type=str,
        nargs="+",
        help="Output file path(s) for reads discarded due to being too short after trimming. Must match number of input files.",
    )
    # discard inline barcode untrimmed reads
    parser.add_argument(
        "-u",
        "--untrimmed-file",
        type=str,
        nargs="+",
        help="Output file path(s) for reads discarded because expected inline barcodes were not found. Must match number of input files.",
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
        "--with-rname-suffix",
        action="store_true",
        help="Indicate if read names have MGI-style suffixes like '/1', '/2', '.1', or '.2' to be stripped.",
    )

    parser.add_argument(
        "--ensure-inline-barcode",
        action="store_true",
        help="If set, reads without the specified inline barcode(s) will be written to the untrimmed files. Requires adapter scheme to have inline barcodes.",
    )

    parser.add_argument(
        "--trim-polyA", action="store_true", help="Enable trimming of polyA/T tails."
    )
    parser.add_argument(
        "--trim-polyA-wo-direction",
        action="store_true",
        help="Trim polyA/T tails regardless of strand information. If not set, trimming depends on strand from adapter scheme.",
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
        "--auto-rc",
        action="store_true",
        help="Automatically reverse complement reads if the library strand (from adapter scheme) is '-'. For paired-end, R1 and R2 will be swapped.",
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

    args = parser.parse_args()

    if args.list_adapters:
        print_builtin_adapters()
        sys.exit(0)

    # Check if input file is provided
    if args.input_file is None:
        logging.error("Input file is required.")
        sys.exit(1)
    elif len(args.input_file) > 2:
        logging.error("Input file can not be more than two.")
        sys.exit(1)

    if args.adapter_name is not None:
        if args.adapter_scheme is not None:
            logging.info("Adapter scheme is provided, ignoring adapter name.")
        else:
            args.adapter_scheme = BUILDIN_ADAPTERS.get(args.adapter_name.upper())
            if args.adapter_scheme is None:
                logging.error(
                    f"Adapter name '{args.adapter_name} not found in built-in adapters."
                )
                # sys.exit(1) # Consider exiting if a bad name is given and no scheme
                # For now, allow fallback to adapter_scheme = adapter_name if not found
                args.adapter_scheme = args.adapter_name
    elif args.adapter_scheme is None:
        logging.error("Adapter scheme or name is required. Use -a or -A.")
        sys.exit(1)
    args.adapter_scheme = args.adapter_scheme.replace(" ", "").upper()

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
    args.short_file = validate_output_file(
        args.short_file, args.input_file, args.output_prefix, "short"
    )

    def _check_with_inline_barcode(s):
        # inline barcode is in bracket () and length > 0
        return re.match(r".*\([ATGCatgc]+\).*", s) is not None

    # Only generate untrimmed file paths if explicitly requested OR if ensure_inline_barcode is true AND the scheme has inline barcodes
    if args.untrimmed_file or (
        args.ensure_inline_barcode and _check_with_inline_barcode(args.adapter_scheme)
    ):
        args.untrimmed_file = validate_output_file(
            args.untrimmed_file, args.input_file, args.output_prefix, "untrimmed"
        )
    else:  # Ensure it's a list of Nones matching input file count if not used
        args.untrimmed_file = [None] * len(args.input_file)

    run_cutseq(args)


if __name__ == "__main__":
    main()
