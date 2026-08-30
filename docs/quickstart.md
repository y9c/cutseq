---
title: Quick Start
nav_order: 2
---

# Quick Start

## Installation

```bash
pip install cutseq
```

## Basic Usage

Trim adapters using a built-in scheme:

```bash
cutseq -A TAKARAV3 test_R1.fq.gz test_R2.fq.gz
```

Or specify a custom adapter scheme in the grammar format:

```bash
cutseq -A "ACACGACGCTCTTCCGATCTXXX-XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC" test_R1.fq.gz test_R2.fq.gz
```

![](https://raw.githubusercontent.com/y9c/cutseq/main/docs/explain_library.png)

The scheme string describes the structure of your library, written 5' to 3':
- Uppercase `ACGT...`: adapters (matched & trimmed)
- Lowercase `acgt...`: inline barcodes (matched & captured)
- `NNNNN` / `N8`: UMI sequences (captured into the read name)
- `XXXXXX` / `X8`: masked sequences (trimmed, not captured)
- `+` `-` `:` in the middle: library strand / R1-R2 split (+ sense, - antisense, : unstranded)

Numeric shorthand (`N8`, `X6`) is equivalent to the expanded run form (`NNNNNNNN`, `XXXXXX`).

## Customizing read names

Captured UMIs/barcodes are appended to the read name by default (legacy naming,
e.g. `@READID_CTATTAAAAA`). Customize with:

```bash
# Custom template (cutadapt brace syntax)
cutseq -A ECLIP10 --name-format '{id}|{cut_prefix}{cut_suffix}' test_R1.fq.gz

# Separator between multiple captured UMIs/barcodes
cutseq -A INLINE --capture-separator ':' test_R1.fq.gz test_R2.fq.gz
```

Defaults reproduce legacy output byte-for-byte.

## Inline barcode auto-detection

Inline barcodes are written in lowercase (`acgt...`). If you write one in
uppercase, it would merge with the adjacent sequencing primer. CutSeq checks
the two outermost adapters against a curated database of Illumina / BGI (MGI)
sequencing primers (`cutseq --list-primers`) and reclassifies any fixed
uppercase sequence adjacent to a recognized primer as an inline barcode, with
a warning. Custom schemes without known primers are never altered.

```bash
# List the known sequencing primers
cutseq --list-primers

# Disable auto-detection
cutseq --no-auto-inline -A "..." test_R1.fq.gz
```

## Discarding reads with a reason tag

A single discard output captures all reads that fail QC, with the reason in
the read name:

```bash
cutseq -A SMALLRNA -d discard_R1.fq.gz discard_R2.fq.gz \
    --min-length 20 --max-n 5 --min-avg-quality 20 \
    test_R1.fq.gz test_R2.fq.gz
```

Discarded reads carry a `reason=...` tag:

- `reason=too_short` — shorter than `--min-length` after trimming
- `reason=too_many_n` — exceeds `--max-n`
- `reason=low_quality` — mean Phred quality below `--min-avg-quality`
- `reason=no_barcode` — missing an expected inline barcode (only with
  `--ensure-inline-barcode`)

`--ensure-inline-barcode` requires the scheme to declare inline barcodes.
When no `-d`/`-O` is given, discard files are auto-named from the input files.

## How it works

A scheme is compiled once into a lazy graph of native cutadapt modifiers and
executed in a single pass. The token-kind registry and phase order are
data-driven, so the engine is easy to extend with new scheme elements. See the [documentation](https://cutseq.yech.science) for details.
