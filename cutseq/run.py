#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-04-19 18:57

import argparse
import copy
import logging
import re
import subprocess
import sys

from cutadapt.adapters import (
    BackAdapter,
    NonInternalBackAdapter,
    NonInternalFrontAdapter,
    PrefixAdapter,
    RightmostFrontAdapter,
)
from cutadapt.files import InputPaths, OutputFiles
from cutadapt.modifiers import (
    AdapterCutter,
    PairedEndRenamer,
    QualityTrimmer,
    SuffixRemover,
    UnconditionalCutter,
)
from cutadapt.pipeline import PairedEndPipeline
from cutadapt.predicates import IsUntrimmed, TooShort
from cutadapt.runners import make_runner
from cutadapt.steps import (
    InfoFileWriter,
    PairedEndFilter,
    PairedEndSink,
    PairedSingleEndStep,
    SingleEndFilter,
    SingleEndSink,
)
from cutadapt.utils import Progress

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
        self.discarded_untrimmed = False
        self.trim_polyA = False
        self.min_length = 20
        self.min_quality = 20
        self.dry_run = False
        self.threads = 1


def run_cutadapt_PE(
    input1, input2, output1, output2, discard1, discard2, barcode, settings
):
    cutadapt = f"cutadapt -j {settings.threads}"
    steps = []
    # step 1: remove adapter on the 5' end, artifact of template switching
    if settings.rname_suffix:
        config_rname = " --strip-suffix '/1' --strip-suffix '/2' --strip-suffix '.1' --strip-suffix '.2'"
    else:
        config_rname = ""
    steps.append(
        f"{cutadapt}{config_rname} -e 0.25 -n 2 -O 10 -g '{barcode.p5.fw};rightmost' -G '{barcode.p7.rc};rightmost' --interleaved {input1} {input2}"
    )
    # step 2: remove adapter on the 3' end, read though in the sequencing
    steps.append(
        f"{cutadapt} -e 0.2 -n 2 -O 3 -a '{barcode.p7.fw}' -A '{barcode.p5.rc}' --interleaved -"
    )
    # step 3: trim inline barcode
    config_inline_args = []
    if barcode.inline5.len > 0:
        config_inline_args.append(f"-g ^{barcode.inline5.fw} -U -{barcode.inline5.len}")
    if barcode.inline3.len > 0:
        config_inline_args.append(f"-G ^{barcode.inline3.rc} -u -{barcode.inline3.len}")
    if barcode.inline5.len + barcode.inline3.len > 0:
        if settings.discarded_untrimmed:
            config_inline_args.append(
                f"--untrimmed-output={discard1} --untrimmed-paired-output={discard2}"
            )
        config_inline = " ".join(config_inline_args)
        steps.append(f"{cutadapt} {config_inline} --interleaved -")
    # step 4: extract UMI
    if barcode.umi5.len + barcode.umi3.len > 0:
        steps.append(
            f"{cutadapt} -u {barcode.umi5.len} -u -{barcode.umi3.len} -U {barcode.umi3.len} -U -{barcode.umi5.len} --rename='{{id}}_{{r1.cut_prefix}}{{r2.cut_prefix}}' --interleaved -"
        )
    else:
        steps.append(f"{cutadapt} --rename='{{id}}' --interleaved -")
    # step 5: mask tail in the RNA, which might be artifact of RT
    if barcode.mask5.len + barcode.mask3.len > 0:
        steps.append(
            f"{cutadapt} -u {barcode.mask5.len} -u -{barcode.mask3.len} -U {barcode.mask3.len} -U -{barcode.mask5.len} --interleaved -"
        )
    # step 6: trim polyA
    if settings.trim_polyA:
        if barcode.strand == "+":
            steps.append(
                f"{cutadapt} -O 6 -e 0.15 -a 'A{{100}}' -G 'T{{100}}' --interleaved -"
            )
        elif barcode.strand == "-":
            steps.append(
                f"{cutadapt} -O 6 -e 0.15 -g 'T{{100}}' -A 'A{{100}}' --interleaved -"
            )
        else:
            logging.info("No strand information provided, skip polyA trimming.")
    # step 7: quality control, remove short reads
    steps.append(
        f"{cutadapt} -q {settings.min_quality} --max-n=0 -m {settings.min_length} --too-short-output={discard1} --too-short-paired-output={discard2} -o {output1} -p {output2} --interleaved -"
    )
    return steps


def run_cutadapt_SE(input1, output1, discard1, barcode, settings):
    cutadapt = f"cutadapt -j {settings.threads}"
    steps = []
    # step 1: remove adapter on the 5' end, artifact of template switching
    if settings.rname_suffix:
        config_rname = " --strip-suffix '/1' --strip-suffix '.1'"
    else:
        config_rname = ""
    steps.append(
        f"{cutadapt}{config_rname} -e 0.25 -n 2 -O 10 -g '{barcode.p5.fw};rightmost' {input1}"
    )
    # step 2: remove adapter on the 3' end, read though in the sequencing
    steps.append(f"{cutadapt} -e 0.2 -n 2 -O 3 -a '{barcode.p7.fw}' -")
    # step 3: trim inline barcode
    config_inline_args = []
    if barcode.inline3.len > 0:
        config_inline_args.append(f"-a {barcode.inline3.fw}$")
    if barcode.inline5.len > 0:
        config_inline_args.append(f"-g ^{barcode.inline5.fw}")
    if barcode.inline5.len + barcode.inline3.len > 0:
        if settings.discarded_untrimmed:
            config_inline_args.append(f"--untrimmed-output={discard1}")

        config_inline = " ".join(config_inline_args)
        steps.append(f"{cutadapt} {config_inline} -")
    # step 4: extract UMI
    if barcode.umi5.len + barcode.umi3.len > 0:
        steps.append(
            f"{cutadapt} -u {barcode.umi5.len} -u -{barcode.umi3.len} --rename='{{id}}_{{cut_prefix}}{{cut_suffix}}' -"
        )
    # step 5: mask tail in the RNA, which might be artifact of RT
    steps.append(f"{cutadapt} -u {barcode.mask5.len} -u -{barcode.mask3.len} -")
    # step 6: trim polyA
    if settings.trim_polyA:
        if barcode.strand == "+":
            steps.append(f"{cutadapt} -O 6 -e 0.15 -a 'A{{100}}' -")
        elif barcode.strand == "-":
            steps.append(f"{cutadapt} -O 6 -e 0.15 -g 'T{{100}}' -")
        else:
            logging.info("No strand information provided, skip polyA trimming.")
    # step 7: quality control, remove short reads
    steps.append(
        f"{cutadapt} -q {settings.min_quality} --max-n=0 -m {settings.min_length} --too-short-output={discard1} -o {output1} -"
    )
    return steps


def run_steps(steps, dry_run=False):
    if dry_run:
        print(
            " |\\\n".join([("    " + s if i > 0 else s) for i, s in enumerate(steps)])
        )
        process = subprocess.run("true", shell=True, capture_output=True)
    else:
        cmd = " | ".join(steps)
        process = subprocess.run(cmd, shell=True, capture_output=True)
    return process.stdout.decode(), process.stderr.decode()


def pipeline_single(input1, output1, short1, barcode, settings):
    trimmer = QualityTrimmer(cutoff_front=0, cutoff_back=15)
    adapter = BackAdapter(sequence="GATCGGAAGA", max_errors=1, min_overlap=3)
    modifiers = [UnconditionalCutter(5), trimmer, AdapterCutter([adapter])]

    inpaths = InputPaths(input1)

    with make_runner(inpaths, cores=settings.threads) as runner:
        outfiles = OutputFiles(
            proxied=settings.threads > 1,
            qualities=runner.input_file_format().has_qualities(),
            interleaved=False,
        )
        steps = [
            # --info-file=info.txt
            PairedSingleEndStep(InfoFileWriter(outfiles.open_text("info.txt"))),
            # -m 10:0
            SingleEndFilter(
                TooShort(settings.min_length),
                outfiles.open_record_writer(short1, interleaved=False),
            ),
            # --discard-untrimmed
            SingleEndFilter(IsUntrimmed()),
            # -o ... -p ...
            SingleEndSink(outfiles.open_record_writer(output1)),
        ]
        pipeline = PairedEndPipeline(modifiers, steps)
        stats = runner.run(pipeline, Progress(), outfiles)
    _ = stats.as_json()
    outfiles.close()


def pipeline_paired(
    input1, input2, output1, output2, short1, short2, barcode, settings
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
        modifiers.append(
            (
                AdapterCutter(
                    [PrefixAdapter(sequence=barcode.inline5.fw, max_errors=0.2)],
                    times=1,
                ),
                UnconditionalCutter(-barcode.inline5.len),
            )
        )
    if barcode.inline3.len > 0:
        modifiers.append(
            (
                UnconditionalCutter(-barcode.inline3.len),
                AdapterCutter(
                    [BackAdapter(sequence=barcode.inline3.rc, max_errors=0.2)],
                    times=1,
                ),
            )
        )
    # TODO: check inline barcode is trimmed in the `info.matches` tag

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
            QualityTrimmer(cutoff_front=0, cutoff_back=15),
            QualityTrimmer(cutoff_front=0, cutoff_back=15),
        )
    )

    inpaths = InputPaths(input1, input2)

    with make_runner(inpaths, cores=settings.threads) as runner:
        outfiles = OutputFiles(
            proxied=settings.threads > 1,
            qualities=runner.input_file_format().has_qualities(),
            interleaved=False,
        )
        steps = [
            # --info-file=info.txt
            # PairedSingleEndStep(InfoFileWriter(outfiles.open_text("info.txt"))),
            # -m 10:10
            PairedEndFilter(
                TooShort(settings.min_length),
                TooShort(settings.min_length),
                outfiles.open_record_writer(discard1, discard2, interleaved=False),
            ),
            # TODO: --max-n=0 support
            # --discard-untrimmed
            # PairedEndFilter( IsUntrimmed(), IsUntrimmed(), pair_filter_mode="any"),
            # -o ... -p ...
            PairedEndSink(outfiles.open_record_writer(output1, output2)),
        ]
        pipeline = PairedEndPipeline(modifiers, steps)
        _stats = runner.run(pipeline, Progress(), outfiles)
        # _ = stats.as_json()
    outfiles.close()


def run_cutseq(args):
    barcode_config = BarcodeConfig(args.adapter_scheme.upper())
    settings = CutadaptConfig()
    if args.with_rname_suffix:
        settings.rname_suffix = True
    if args.discarded_untrimmed:
        settings.discarded_untrimmed = True
    if args.trim_polyA:
        settings.trim_polyA = True
    settings.threads = args.threads
    settings.min_length = args.min_length
    settings.dry_run = args.dry_run
    if len(args.input_file) == 1:
        pipeline_single(
            args.input_file[0],
            args.output_file[0],
            args.discard_file[0],
            barcode_config,
            settings,
        )
    else:
        pipeline_paired(
            args.input_file[0],
            args.input_file[1],
            args.output_file[0],
            args.output_file[1],
            args.discard_file[0],
            args.discard_file[1],
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
    # output file can be number of files matching the input files, if not provided it will generate based on the output suffix,
    # if no output suffix provided it will generate based on the input file name
    parser.add_argument(
        "-a",
        "--adapter-scheme",
        type=str,
        help="Adapter sequence configuration.",
    )
    parser.add_argument("-A", "--adapter-name", type=str, help="Built-in adapter name.")
    parser.add_argument(
        "-O",
        "--output-suffix",
        type=str,
        help="Output file suffix for keep trimmed data.",
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
        "-S",
        "--short-file",
        type=str,
        nargs="+",
        help="Output file path for discarded too short data.",
    )
    parser.add_argument(
        "-m",
        "--min-length",
        type=int,
        default=20,
        help="Minimum length of the reads to keep.",
    )
    parser.add_argument(
        "-q",
        "--min-quality",
        type=int,
        default=20,
        help="Minimum quality of the read tails in the reads to keep..",
    )

    parser.add_argument(
        "--with-rname-suffix",
        action="store_true",
        help="R1 and R2 suffix cotains suffix. MGI platform.",
    )
    parser.add_argument(
        "-U",
        "--untrimmed-file",
        type=str,
        nargs="+",
        help="Output file path for discarded reads without inline barcode.",
    )
    parser.add_argument("--trim-polyA", action="store_true", help="Trim polyA tail.")

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
    args = parser.parse_args()

    if args.adapter_name is not None:
        if args.adapter_scheme is not None:
            logging.info("Adapter scheme is provided, ignore adapter name.")
        else:
            if args.adapter_name.upper() == "TAKARAV2":
                args.adapter_scheme = "ACACGACGCTCTTCCGATCTX<XXXAGATCGGAAGAGCACACGTC"
            elif args.adapter_name.upper() == "STRANDED":
                args.adapter_scheme = "ACACGACGCTCTTCCGATCTX<XXXAGATCGGAAGAGCACACGTC"
            elif args.adapter_name.upper() == "TAKARAV3":
                args.adapter_scheme = (
                    "ACACGACGCTCTTCCGATCTXXX<XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"
                )
            else:
                logging.error("Adapter name is not valid.")
                sys.exit(1)
    elif args.adapter_scheme is None:
        logging.error("Adapter scheme or name is required.")
        sys.exit(1)

    if len(args.input_file) > 2:
        raise ValueError("Input file can not be more than two.")

    if args.output_file:
        if len(args.output_file) != len(args.input_file):
            raise ValueError("Output file should be same as input file.")
    elif args.output_suffix:
        if len(args.input_file) == 1:
            args.output_file = [args.output_suffix + "_trimmed_R1.fastq.gz"]
        else:
            args.output_file = [
                args.output_suffix + "_trimmed_R1.fastq.gz",
                args.output_suffix + "_trimmed_R2.fastq.gz",
            ]
    else:
        if len(args.input_file) == 1:
            args.output_file = [
                remove_fq_suffix(args.input_file[0]) + "_trimmed_R1.fastq.gz",
            ]
        else:
            args.output_file = [
                remove_fq_suffix(args.input_file[0]) + "_trimmed_R1.fastq.gz",
                remove_fq_suffix(args.input_file[1]) + "_trimmed_R2.fastq.gz",
            ]

    if args.short_file:
        if len(args.short_file) != len(args.input_file):
            raise ValueError("short file should be same as input file.")
    elif args.output_suffix:
        if len(args.input_file) == 1:
            args.short_file = [args.output_suffix + "_shorted_R1.fastq.gz"]
        else:
            args.short_file = [
                args.output_suffix + "_shorted_R1.fastq.gz",
                args.output_suffix + "_shorted_R2.fastq.gz",
            ]
    else:
        if len(args.input_file) == 1:
            args.short_file = [
                remove_fq_suffix(args.input_file[0]) + "_shorted_R1.fastq.gz",
            ]
        else:
            args.short_file = [
                remove_fq_suffix(args.input_file[0]) + "_shorted_R1.fastq.gz",
                remove_fq_suffix(args.input_file[1]) + "_shorted_R2.fastq.gz",
            ]
    run_cutseq(args)


if __name__ == "__main__":
    main()
