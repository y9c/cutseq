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

The customized scheme can be explained by the diagram above. The scheme is a
**top-strand molecular map**: you write the library as one DNA strand, 5' to 3',
with no spaces:

- The part **left of the insert marker** (`+` `-` `:`) is read by Read 1 as-is
  (R1 sequences the top strand 5' -> 3').
- The part **right of the insert marker** is the molecule's 3' continuation,
  written top-strand 5' -> 3' (from the insert toward p7). Read 2 sequences the
  bottom strand, so the engine matches it as the **reverse complement** of that
  part automatically — e.g. a right-hand `AGATCGGAAGAGCACACGTC` is matched on
  Read 2 as `GACGTGTGCTCTTCCGATCT`.

- The outermost parts on both ends are the Illumina adapters (uppercase `ACGT...`).
- The UMI sequence is the random sequence represented by `N` (e.g. `NNNNNNNN` / `N8`), captured into the read name.
- The masked sequences are represented by `X` (e.g. `XXXXXX` / `X6`), trimmed but not captured. These can be random tails from template switching or other artifacts.
- Homopolymer tails are trimmed with the dot form `AAA...AAA` / `TTT...TTT` / `GGG...GGG`. The 5'/3' direction is **auto-detected from the scheme layout**: a run at a read start (first token, or right after the outer sequencing adapter/primer) trims the 5′ (leftmost-anchored — e.g. a template-switching poly‑G stretch), while a run elsewhere trims the 3′ (e.g. a poly‑A tail). `B10` style (`G10`) sets the minimum run length.
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
customize this with `--rename` (`--name-format` is an accepted alias). It
supports cutadapt's brace variables (`{id}`, `{header}`, `{comment}`,
`{cut_prefix}`, `{cut_suffix}`, `{adapter_name}`, `{match_sequence}`, `{rc}`)
plus **positional captures** — the individual UMIs and inline barcodes in
scheme order:

```bash
# Label each captured part; {1},{2},{3} are the 1st/2nd/3rd capture
cutseq -A INLINE --rename '{id}_BC1:{1}_BC2:{2}_umi:rc({3})' in_R1.fq.gz in_R2.fq.gz
```

Transform **functions** can wrap any capture and nest, and are case-insensitive:

- `rc(x)` reverse complement (also `rev(x)` reverse, `comp(x)` complement only)
- `upper(x)`, `lower(x)`, `len(x)`
- `canon(x)` — canonical UMI form `min(x, rc(x))` for UMI collapsing
- `left(x,k)`, `right(x,k)`, `slice(x,a,b)` — substring helpers
  (`slice` is 1-based, inclusive)
- examples: `upper(RC({1}))`, `rc(upper({2}))`, `left({2},6)`, `slice({1},1,4)`

In paired mode `{r1.1}` / `{r2.1}` force a specific read; an unprefixed capture
resolves to its *anchor* read (left-side captures → R1, right-side → R2), so
both mates carry the same extracted value.

Use `--rename` to reproduce the old `--capture-separator ':'` behavior too:
`--rename '{id}:{1}:{2}'` inserts `:` between the captured parts.

When no `--rename` is given, the default naming (captures appended with `_`)
reproduces legacy output byte-for-byte, so existing pipelines do not change.

## Complex example: spatial barcode arm (DBiT-seq)

Spatial libraries put several fixed-length barcodes on one arm. In **DBiT-seq**
(deterministic barcoding in tissue, 50×50 grid), the barcode read (R2) walks:

```
handle(22) | BarcodeB(8) | linker(30) | BarcodeA(8) | linker(30) | UMI(12) | polyT | cDNA
```

while Read 1 is the cDNA read, preceded by a template-switch oligo (TSO).
Everything fits in one command:

```bash
cutseq -A DBITSEQ \
       -R '{id}_BCB:{3}_BCA:{2}_UMI:{1}' \
       -O mysample R1.fq.gz R2.fq.gz
```

Built-in scheme (`cutseq -A DBITSEQ`):

```
AAGCAGTGGTATCAACGCAGAGT : N12 AGTCGTACGCCGATGCGAAACATCGGCCAC
N8 CGAATGCTCTGGCCTCTCAAGCACGTGGAT N8 AGATGCGAGAAGCCAACGCTTG
```

- The `TSO` (left of `:`) is R1's real read-5′ scaffold — R1 actually starts
  with it, downstream of the R1 sequencing primer — so it is trimmed off R1.
- The barcode arm (right of `:`) is walked on R2 as its reverse complement:
  handle trimmed, then `BarcodeB` (8 nt, 1st capture → `{3}`), linker 2 trimmed,
  `BarcodeA` (8 nt, 2nd capture → `{2}`), linker 1 trimmed, `UMI` (12 nt, 3rd
  capture → `{1}`). Each captured part is written into the read header as
  `@READID_BCB:…_BCA:…_UMI:…`.
- `M6AARTR` is kept as a deprecated alias; `--r1-primer` / `--r2-primer` are
  **informational only** — reads start downstream of the priming site, so the
  sequencing primers themselves are never trimmed.

The same library as a fully inline custom scheme:

```bash
cutseq -A "AAGCAGTGGTATCAACGCAGAGT:N12AGTCGTACGCCGATGCGAAACATCGGCCACN8CGAATGCTCTGGCCTCTCAAGCACGTGGATN8AGATGCGAGAAGCCAACGCTTG" \
       -R '{id}_BCB:{3}_BCA:{2}_UMI:{1}' -O mysample R1.fq.gz R2.fq.gz
```

Full construct, top strand 5′→3′:

```
P5(29) | i5(8) | R1-primer(34) | TSO(22) | insert | polyA/T | UMI(12)
| rc(L1)(30) | BarcodeA(8) | rc(L2)(30) | BarcodeB(8) | rc(handle)(22)
| rc(R2-primer)(34) | i7(8) | rc(P7)(24)
```

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
