---
title: Spatial (DBiT-seq) libraries
nav_order: 5
---

# Parsing spatial (DBiT-seq) libraries

This page is the introduction to handling **spatial sequencing** data with
CutSeq: what the library construct looks like, which read carries what, the
one-line command that trims everything and captures the spatial barcodes into
the read names, and how to interpret the QC that results from it.

It is written for the **DBiT-seq** scheme (`cutseq -A DBITSEQ`, spatial RNA,
50&times;50 barcode grid), but the same walk-through applies to any
barcode-arm library expressed as a CutSeq scheme.

---

## 1. The library construct

DBiT-seq puts a row of fixed-length oligos — the **barcode arm** — on one side
of each molecule. The full construct, top strand 5&prime;&rarr;3&prime;:

```
P5(29) | i5(8) | R1-primer(34) | TSO(23) | insert | polyA/T | UMI(12)
       | rc(L1)(30) | BarcodeA(8) | rc(L2)(30) | BarcodeB(8)
       | rc(handle)(22) | rc(R2-primer)(34) | i7(8) | rc(P7)(24)
```

The sequencing primers just upstream of each read are **never part of the
reads** — sequences start downstream of the priming site. So:

| Read | What it actually sequences |
|------|----------------------------|
| **R1** (from the P5 side) | **TSO** (23 nt template-switch oligo) + cDNA **insert** (variable length), often with a 5&prime; G-stretch right after the TSO |
| **R2** (from the P7 side) | **barcode arm** then insert read-through: handle(22) → **BarcodeB**(8) → linker2(30) → **BarcodeA**(8) → linker1(30) → **UMI**(12) → polyA/T → 3&prime; end of the insert |

{: .note }
The insert is on R1 (5&prime;&rarr;3&prime;, top strand). R2 only sees the
**3&prime; tail** of the insert that reads through past the arm — and it is
capped by the sequencer read length (see
[§5 – interpreting read lengths](#5-interpreting-the-read-length-distribution)).

## 2. One command to parse everything

```bash
cutseq -A DBITSEQ \
  -R '{id}_BCB:{3}_BCA:{2}_UMI:{1}' \
  -O mysample \
  mysample_R1.fastq.gz mysample_R2.fastq.gz
```

Fast, single-threaded by default (`-t N` to parallelize). That one line:

* trims the **TSO + 5&prime; G-stretch off R1**,
* walks the **barcode arm on R2** (handle, both linkers, poly tail),
* **captures BarcodeB, BarcodeA and the UMI** and writes them into every read
  name, e.g.:

```
@LH00699:99:2533WTLT4:6:1101:36114:1055_BCB:AGAGTCAA_BCA:ACCTCCAA_UMI:TGCACACCAACC
```

### What the built-in `DBITSEQ` scheme means

```
AAGCAGTGGTATCAACGCAGAGTGAATGGG...GGG : N12 AGTCGTACGCCGATGCGAAACATCGGCCAC
                                            N8  CGAATGCTCTGGCCTCTCAAGCACGTGGAT
                                            N8  AGATGCGAGAAGCCAACGCTTG
```

| Token | Side | Meaning in CutSeq |
|-------|------|-------------------|
| `AAGCAGTGGTATCAACGCAGAGT` | left (`:`⇒R1, as-is) | TSO — R1&rsquo;s real read-5&prime; scaffold, trimmed off |
| `GAATGGG...GGG` | R1 | auto-detected **5&prime; poly-Run** (G-stretch / GAATGGG) leftmost-anchored — trims the template-switch remnant; `G<k>` tunes the minimum run |
| `N12` | right (matched on R2) | **UMI** capture &rarr; `{1}` |
| first 30-nt run | right | linker 1, trimmed |
| `N8` | right | **BarcodeA** capture &rarr; `{2}` |
| second 30-nt run | right | linker 2, trimmed |
| `N8` | right | **BarcodeB** capture &rarr; `{3}` |
| final 22-nt run | right | handle (rc form), trimmed |

The `-R '{id}_BCB:{3}_BCA:{2}_UMI:{1}'` template is read-name renaming —
see the [rename page](rename.md) for the full syntax (transforms, `rc()`, per-side
captures, etc.).

## 3. Outputs

With `-O mysample` you get (input file names are used if `-O` is omitted):

```
mysample_trimmed_R1.fastq.gz   # insert read, barcodes in header, TSO/G-stretch removed
mysample_trimmed_R2.fastq.gz   # barcode arm trimmed, barcodes in header, UMI in header
mysample_discard_R1.fastq.gz   # rejected reads (both mates kept together)
mysample_discard_R2.fastq.gz
```

Rejected reads carry a `reason=` tag in the name: `too_short` (below
`-m/--min-length`, default 20), `too_many_n` (`--max-n`), `low_quality`
(`--min-avg-quality`), `no_barcode` (only with `--ensure-inline-barcode`).
Add `--json-file mysample.json` for a full trimming report with
counts/statistics per step.

{: .note }
Because trimming is performed on paired, in-order reads, R1 (insert) and R2
(barcode) mates stay in sync — the barcode in an R2 header belongs to the
insert sequence in the same R1 record.

## 4. Going from trimmed reads to a spatial expression matrix

1. **Parse the header** — split each read name on `_BCB:`, `_BCA:`, `_UMI:`
   and read the 8+8+12 nt.
2. **Map (BarcodeB, BarcodeA) → pixel** with the design&rsquo;s barcode whitelist
   (the two 8-nt sets used to build the 50&times;50 grid; allow 1–2 nt
   Levenshtein distance if the protocol tolerates it).
3. **Align R1** to the reference and assign transcripts/features per pixel
   (m6A-ARTR-style methods typically report per-feature, per-pixel read counts).
4. **Deduplicate by UMI** — collapse reads sharing the same pixel + feature
   (+ alignment position) and the same 12-nt UMI.

A minimal header parser / demultiplexer (paired trimmed R1+R2, barcodes in the
R1 name):

```python
import re, gzip

PAT = re.compile(r"_BCB:([ACGT]{8})_BCA:([ACGT]{8})_UMI:([ACGT]{12})")

def records(fq):
    with gzip.open(fq, "rt") as f:
        lines = [ln.strip() for ln in f]
    for i in range(0, len(lines), 4):
        yield lines[i], lines[i + 1]                # (name, sequence)

with open("pixels.tsv", "w") as out:
    for (r1n, r1s), (r2n, r2s) in zip(
            records("mysample_trimmed_R1.fastq.gz"),
            records("mysample_trimmed_R2.fastq.gz")):
        m = PAT.search(r1n)
        if not m:
            continue
        bcb, bca, umi = m.groups()
        out.write(f"{r1n.split()[0]}\t{bcb}\t{bca}\t{umi}\t{r1s}\n")
```

## 5. Interpreting the read-length distribution (the R2 “41-nt peak”)

After trimming, R1 lengths reflect the real insert-length spread (here ~9–123 nt,
peaking at ~11–15). **R2 is different**: the fixed barcode arm occupies 110 nt
(handle 22 + BCB 8 + linker2 30 + BCA 8 + linker1 30 + UMI 12) of the sequenced
read, so the arm + UMI + read-through has to fit in the read.

For a 151-cycle run:

```
151 (read length) − 110 (arm: handle22 + BCB8 + L2 30 + BCA8 + L1 30 + UMI12) = 41 nt
```

Every molecule whose insert is **longer than ~41 nt on the R2 side** therefore
hits the physical end of the read and piles up at exactly **41 nt** after
trimming (`raw[-41:]` in 99% of cases). In practice this can be ~18% of all R2
reads — it is a **read-truncation artefact, not a biological 41-nt fragment**:
the sequences are fully diverse and their R1 mates are long.

{: .warning }
Do **not** read R2&rsquo;s length axis as fragment length. R2 is capped by the
sequencer, not by the insert. If you need the full 3&prime; R2-side sequence
read length 151→250 (`--r2-primer` is informational only), or rely on **R1**
for insert biology.

Correspondingly, in a per-position **base-composition &times; depth** pileup:
* **R1** — a clean, slightly composition-skewed stack that runs out around the
  longest insert lengths;
* **R2** — strong constant composition across the arm-bearing positions, then a
  tail of read-through cDNA; expect the depth to collapse once reads that
  reached the 151-nt end stop contributing.

The built-in defaults are gentle (`-m 20`, `-q 20`, no `--max-n` /
`--min-avg-quality` by default); set them from your QC thresholds and check the
`reason=` distribution in the discard files.

## 6. Notes & extras

* Legacy name: `-A M6AARTR` still works as an alias for `-A DBITSEQ` — prefer
  `-A DBITSEQ`.
* The same library as a fully inline custom scheme (no built-in needed):
  `cutseq -A "AAGCAGTGGTATCAACGCAGAGTGAATGGG...GGG:N12AGTCGTACGCCGATGCGAAACATCGGCCACN8CGAATGCTCTGGCCTCTCAAGCACGTGGATN8AGATGCGAGAAGCCAACGCTTG" -R '{id}_BCB:{3}_BCA:{2}_UMI:{1}' ...`
* Dry-run inspection: `cutseq -A DBITSEQ --graph-vertical R1.fq.gz R2.fq.gz`
  prints the parsed trimming steps without writing any output — useful to
  confirm the arm layout before running the real job.
* `--r1-primer` / `--r2-primer` document the sequencing primers but are never
  trimmed from reads (reads start downstream of them).
* See [Quick Start](quickstart.md), [Adapter schemes](adapters.md) and
  [Read-name renaming](rename.md) for more.
