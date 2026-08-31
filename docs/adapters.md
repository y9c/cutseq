---
title: Adapter Schemes
nav_order: 3
---

# Adapter Schemes

CutSeq supports a variety of built-in adapter schemes for common NGS library types. You can list all available schemes in your terminal with:

```bash
cutseq --list-adapters
```

Use the adapter name (or a custom grammar scheme string) with `-A/--adapter-scheme`. Note `-a` is accepted as an alias of `-A`.

## Example: Built-in Schemes

- **SMALLRNA**: Small RNA libraries, double ligation, forward orientation
- **INLINE**: Custom barcoded libraries, dual UMI, inline barcode
- **TAKARAV3**: SMARTer Stranded Total RNA-Seq Kit v3
- **STRANDED**: Stranded RNA libraries

See below for a comprehensive guide to each supported adapter pattern, including copyable scheme blocks and usage notes.

---

## Inline barcode auto-detection

Inline barcodes are written in lowercase (`acgt...`). If you write one in uppercase, it merges with the adjacent sequencing primer into one adapter run. CutSeq checks the two outermost adapters of a scheme against a curated database of Illumina / BGI (MGI) sequencing primers (`cutseq --list-primers`) and reclassifies any fixed uppercase sequence adjacent to (or between) recognized primers as an inline barcode, with a warning. Custom schemes without known primers are never altered. Disable with `--no-auto-inline`.

## Discarding reads with a reason tag

A single discard output captures all reads that fail QC, with the reason stored in the read name (`reason=too_short`, `reason=too_many_n`, `reason=low_quality` or `reason=no_barcode`):

```bash
cutseq -A SMALLRNA -d discard_R1.fq.gz discard_R2.fq.gz \
    --min-length 20 --max-n 5 --min-avg-quality 20 \
    test_R1.fq.gz test_R2.fq.gz
```

- `reason=too_short` — shorter than `--min-length` after trimming
- `reason=too_many_n` — exceeds `--max-n`
- `reason=low_quality` — mean Phred quality below `--min-avg-quality`
- `reason=no_barcode` — missing an expected inline barcode (only with `--ensure-inline-barcode`, which requires inline barcodes in the scheme)

When no `-d`/`-O` is given, discard files are auto-named from the input files.

---

### SMALLRNA (Small RNA Libraries)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: AGTTCTACAGTCCGACGATC+AGATCGGAAGAGCACACGTC" data-scheme="AGTTCTACAGTCCGACGATC+AGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGTTCTACAGTCCGACGATC</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Used for small RNA sequencing
- Double ligation method
- Forward orientation
- Optional 2nt trimming on both ends for quality

---
### SMRNA (Small RNA, Legacy version)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: AGTTCTACAGTCCGACGATC+TGGAATTCTCGGGTGCCAAG" data-scheme="AGTTCTACAGTCCGACGATC+TGGAATTCTCGGGTGCCAAG"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGTTCTACAGTCCGACGATC</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">TGGAATTCTCGGGTGCCAAG</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Legacy version for small RNA sequencing
- Refer to library preparation protocol for details.

---
### INLINE (Custom Barcoded Libraries)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: AGTTCTACAGTCCGACGATCNNNNN+NNNNNatcacgAGATCGGAAGAGCACACGTC" data-scheme="AGTTCTACAGTCCGACGATCNNNNN+NNNNNatcacgAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGTTCTACAGTCCGACGATC</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNN</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNN</span><span style="background-color: #FFD700; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">atcacg</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Used for libraries with inline barcodes
- Dual UMI design (5nt each)
- Forward orientation
- Contains fixed inline barcode (ATCACG)

---
### TAKARAV2 (Takara V2)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXXX-XXXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXXX-XXXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XXX)
- Refer to Takara V2 library preparation protocol for details.

---
### STRANDED (Stranded RNA)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTX-XXXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTX-XXXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">X</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (X, XXX)
- Refer to stranded RNA library preparation protocol for details.

---
### UNSTRANDED (Unstranded RNA/DNA)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXX:XXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXX:XXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">:</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XX)
- Refer to library preparation protocol for details.

---
### TAKARAV3 (Takara V3)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXXX-XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXXX-XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXX</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XXX, XXXXXX)
- Contains UMI sequences (NNNNNNNN)
- Refer to Takara V3 library preparation protocol for details.

---
### ECLIP6 (eCLIP (6nt UMI))

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXX-XNNNNNNAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXX-XNNNNNNAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">X</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XX, X)
- Contains UMI sequences (NNNNNN)
- Used for eCLIP experiments with a 6-nucleotide UMI.

---
### ECLIP10 (eCLIP (10nt UMI))

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXX-XNNNNNNNNNNAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXX-XNNNNNNNNNNAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">X</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XX, X)
- Contains UMI sequences (NNNNNNNNNN)
- Used for eCLIP experiments with a 10-nucleotide UMI.

---
### SACSEQ (SAC-Seq)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCT-XXXXNNNNNNAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCT-XXXXNNNNNNAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXX</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XXXX)
- Contains UMI sequences (NNNNNN)
- Used for SAC-Seq (Single-cell RNA-Seq) or similar protocols.

---
### SACSEQV3 (SAC-Seq V3)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTNNNNNNNNX+XXNNNNNNNNAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTNNNNNNNNX+XXNNNNNNNNAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNN</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">X</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XX</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (X, XX)
- Contains dual UMI sequences (NNNNNNNN on both sides)
- Refer to SAC-Seq V3 library preparation protocol for details.

---
### XGENRNA (xGen RNA)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXXXXXX-XXXXXXXXXXXXXXXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXXXXXX-XXXXXXXXXXXXXXXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXXXXXXXXXXX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains extensive masked regions (XXXXXX, XXXXXXXXXXXXXXX)
- Refer to xGen RNA Library Preparation Kit protocol for details.

---
### ILLUMINARNA (Illumina Stranded RNA (Nextera-based))

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: AGATGTGTATAAGAGACAG-CTGTCTCTTATACACATCT" data-scheme="AGATGTGTATAAGAGACAG-CTGTCTCTTATACACATCT"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATGTGTATAAGAGACAG</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">CTGTCTCTTATACACATCT</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Uses Nextera-style adapter sequences.

---
### DSLIGATION (dsDNA Ligation)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCT+AGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCT+AGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Basic dsDNA ligation with A-tailing
- No UMIs or special trimming needed

---
### XGENMETHY (xGen Methyl-Seq)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXX+XXXXXXXXXXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXX+XXXXXXXXXXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXXXXXX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XX, XXXXXXXXXX)
- Designed for xGen Methyl-Seq library kits.

---
### XGENSNMC (xGen snmC-Seq)

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXXXXXX+XXXXXXXXXXXXXXXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXXXXXX+XXXXXXXXXXXXXXXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXXXXXXXXXXX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains extensive masked regions (XXXXXX, XXXXXXXXXXXXXXX)
- Designed for xGen single-nucleus methylC-Seq (snmC-Seq) library kits.

---
### PBAT (PBAT (Post-Bisulfite Adapter Tagging))

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: ACACGACGCTCTTCCGATCTXXXXXX-XXXXXXAGATCGGAAGAGCACACGTC" data-scheme="ACACGACGCTCTTCCGATCTXXXXXX-XXXXXXAGATCGGAAGAGCACACGTC"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">ACACGACGCTCTTCCGATCT</span><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXX</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">-</span></div><span style="background-color: #DCDCDC; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">XXXXXX</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATCGGAAGAGCACACGTC</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Contains masked regions (XXXXXX on both sides)
- Used for PBAT-style bisulfite sequencing libraries.

---
### NEXTERA (Nextera (General))

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: AGATGTGTATAAGAGACAG+CTGTCTCTTATACACATCT" data-scheme="AGATGTGTATAAGAGACAG+CTGTCTCTTATACACATCT"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATGTGTATAAGAGACAG</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">+</span></div><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">CTGTCTCTTATACACATCT</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Standard Nextera transposase-based library preparation.
- No UMIs by default in this basic scheme.
- Usually used for ATAC-seq libraries.

---
### M6AARTR (m6A-ARTR-DBiT (spatial))

<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;"><div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: AAGCAGTGGTATCAACGCAGAGT:N10AGTCGTACGCCGATGCGAAACATCGGCCACN8CGAATGCTCTGGCCTCTCAAGCACGTGGATN8AGATGCGAGAAGCCAACGCTTG" data-scheme="AAGCAGTGGTATCAACGCAGAGT:N10AGTCGTACGCCGATGCGAAACATCGGCCACN8CGAATGCTCTGGCCTCTCAAGCACGTGGATN8AGATGCGAGAAGCCAACGCTTG"><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AAGCAGTGGTATCAACGCAGAGT</span><div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">:</span></div><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGTCGTACGCCGATGCGAAACATCGGCCAC</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">CGAATGCTCTGGCCTCTCAAGCACGTGGAT</span><span style="background-color: #B2EBF2; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">NNNNNNNN</span><span style="background-color: #A8E6CF; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">AGATGCGAGAAGCCAACGCTTG</span></div><div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div></div>
- Spatial RNA library (m6A-ARTR / DBiT).
- R1: TSO + cDNA; the barcode arm is on R2 (written in rc form after ':'),
- per the molecule's 5'->3' order.
- R1's read 5' scaffold is the TSO (template-switch oligo; R1's first
- sequenced bases, downstream of the R1-primer site) and is what gets
- trimmed by the scheme's left token (S0 R1).
- The R1/R2 sequencing primers (TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG and
- GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG) anneal UPSTREAM of each read; reads
- start downstream of them, so the primers are never part of the reads and
- are not trimmed. --r1-primer / --r2-primer are informational only.
- Full construct (5'->3'): P5 | i5(8) | R1-primer(34) | TSO(22) | insert |
- polyA/T | UMI | rc(L1) | barcodeA(8) | rc(L2) | barcodeB(8) | rc(handle) |
- rc(R2-primer) | i7(8) | rc(P7).

---
<script>(function() {  function showTooltip(el) {    var tooltip = el.parentElement.querySelector(".scheme-raw-tooltip");    if (tooltip) {      tooltip.style.display = "block";      setTimeout(function() { tooltip.style.display = "none"; }, 1200);    }  }  document.querySelectorAll(".copy-scheme-raw").forEach(function(block) {    block.addEventListener("mouseenter", function() {      block.style.boxShadow = "0 0 0 2px #FF6F61";    });    block.addEventListener("mouseleave", function() {      block.style.boxShadow = "";    });    block.addEventListener("click", function(e) {      var scheme = block.getAttribute("data-scheme");      if (navigator.clipboard) {        navigator.clipboard.writeText(scheme).then(function() {          showTooltip(block);        });      } else {        var textarea = document.createElement("textarea");        textarea.value = scheme;        document.body.appendChild(textarea);        textarea.select();        document.execCommand("copy");        document.body.removeChild(textarea);        showTooltip(block);      }    });  });})();</script>