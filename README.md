# ✂️ CutSeq

[![Pypi Releases](https://img.shields.io/pypi/v/cutseq.svg)](https://pypi.python.org/pypi/cutseq)
[![Downloads](https://pepy.tech/badge/cutseq)](https://pepy.tech/project/cutseq)

CutSeq is a tool that provides an efficient wrapper for the cutadapt tool, which is powerful in handling various types of NGS libraries.
Due to the complexities involved in NGS library preparation methods, mutiple operations are necessary to process sequencing reads correctly.

Take _SMARTer® Stranded Total RNA-Seq Kit v3_ as an example, at least 9 operations are required.

![](https://raw.githubusercontent.com/y9c/cutseq/main/docs/takaraV3.png)

For **Read 1**:

1.  Remove the Illumina p7 adapter from the end of the sequence.
2.  Remove 14 nt (8+3+3) at the rightmost position of the sequence, representing UMI and linker sequence from the beginning of read 2. This is required when the library insert size is shorter than the sequencing length.
3.  Remove poly-T sequences at the beginning of the sequence (read 1 is oriented in reverse to the RNA, hence a polyA tail appears as a leading polyT sequence).
4.  Remove low-quality bases from right to left.

For **Read 2**:

5.  Remove the reverse complement Illumina p5 adapter from the end of the sequence.
6.  Extract the 8 nt UMI sequence from the beginning of the sequence and append it to the read name for downstream analysis.
7.  Mask a 6 nt linker sequence at the leftmost position immediately after clipping the UMI sequence.
8.  Remove poly-A sequences at the end of the read.
9.  Remove low-quality bases from right to left.

These operations must be performed in the **correct order**. The limitations of the cutadapt tool make it challenging to configure these operations in a single command, often leading to errors unnoticed in some publications.

---

To solve this by using cutadapt, we can run multiple cutadpat insitent sequentially or pipe multiple commands together. But this waste lots of IO and computational resource. I am thinking there a more eligent API to make things easy. Then comes this toy project.
-- **What you need is only one parameter which spcific what the library would looks like.**

CutSeq overcomes these limitations by enabling multiple operations in a automatical manner to ensure accuracy and efficiency.

## How to install?

```bash
pip install cutseq
```

## How to use?

Execute adapter trimming by providing a single parameter and your input files:

```bash
cutseq -A TAKARAV3 test_R1.fq.gz test_R2.fq.gz
```

Alternatively, you can specify a custom adapter scheme in the grammar format:

`cutseq -A "ACACGACGCTCTTCCGATCTXXX-XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"`

![](https://raw.githubusercontent.com/y9c/cutseq/main/docs/explain_library.png)

The customized scheme can be explained by the diagram above. The scheme is
written 5' to 3' (R1 direction) with no spaces required:

- The outermost parts on both ends are the Illumina adapters (uppercase `ACGT...`).
- The UMI sequence is the random sequence represented by `N` (e.g. `NNNNNNNN` / `N8`), captured into the read name.
- The masked sequences are represented by `X` (e.g. `XXXXXX` / `X6`), trimmed but not captured. These can be random tails from template switching or other artifacts.
- Inline barcodes are written in lowercase (`acgt...`) and are matched and captured.
- The center parts are the actual library sequence, split by `+`, `-`, or `:`:
  - `+` means the library is in the forward (sense) orientation,
  - `-` means reverse (antisense) orientation,
  - `:` means the library orientation is unknown (unstranded).
- A 3' read-through adapter (the adapter that appears on the opposite read)
  is trimmed automatically from the other end of each read.

Numeric shorthand (`N8`, `X6`) is equivalent to the expanded run form
(`NNNNNNNN`, `XXXXXX`).

## Customizing read names

By default, captured UMIs/barcodes are appended to the read name with `_`
(e.g. `@READID_CTATTAAAAA`), exactly as the legacy engine named reads. You can
customize this:

```bash
# Custom template using cutadapt's brace syntax
cutseq -A ECLIP10 --name-format '{id}|{cut_prefix}{cut_suffix}' in_R1.fq.gz

# Insert a separator between multiple captured UMIs/barcodes
cutseq -A INLINE --capture-separator ':' in_R1.fq.gz in_R2.fq.gz
```

The default `--name-format` and `--capture-separator` values reproduce legacy
output byte-for-byte, so existing pipelines do not need to change.

## Inline barcodes and auto-detection

Inline barcodes are written in lowercase (`acgt...`) and are matched and
trimmed. If you write a barcode in uppercase by mistake, it would merge with
the adjacent sequencing primer into one adapter run. CutSeq detects this:

- The two outermost adapters of a scheme are checked against a curated
  database of Illumina / BGI (MGI) sequencing primers (see `cutseq --list-primers`).
- Any fixed uppercase sequence adjacent to (or between) recognized primers is
  reclassified as an inline barcode, with a warning.
- Genuinely custom schemes (no known primers at the ends) are never altered.

Disable with `--no-auto-inline`.

## Discarding reads with a reason tag

A single discard output captures all reads that fail QC, with the reason
stored in the read name:

```bash
cutseq -A SMALLRNA -d discarded_R1.fq.gz discarded_R2.fq.gz in_R1.fq.gz in_R2.fq.gz
```

Discarded reads carry a `reason=...` tag in their name:

- `reason=too_short` — shorter than `--min-length` after trimming,
- `reason=too_many_n` — exceeds `--max-n`,
- `reason=low_quality` — mean Phred quality below `--min-avg-quality`,
- `reason=no_barcode` — missing an expected inline barcode (only with
  `--ensure-inline-barcode`).

`--ensure-inline-barcode` routes reads that lack the scheme's inline barcodes
to the discard output; it requires the scheme to declare inline barcodes.
When no `-d`/`-O` is given, discard files are auto-named from the input files.

## How it works

CutSeq compiles a library scheme once into a lazy graph of native cutadapt
modifiers (adapter cutters, renamers, quality trimmers) and runs them in a
single cutadapt pass — no intermediate I/O. The modifier chain is data-driven:
each grammar token kind maps to a 5' / 3' / single-end emitter, so new token
kinds are added with one registry entry. Output from the grammar engine is
verified to match the legacy engine exactly across all built-in schemes in
both paired- and single-end mode.

More details can be found in the [document](https://cutseq.yech.science)

## TODO

[ ] support more library scheme
