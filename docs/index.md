---
title: Home
layout: home
nav_order: 1
description: "CutSeq: automatically cut adapter / barcode / UMI from NGS data."
permalink: /
---

# ✂️ CutSeq
{: .fs-9 }

Automatically cut adapter / barcode / UMI from NGS data, in the right order.
{: .fs-6 .fw-300 }

[Get started](quickstart.md){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Adapter schemes](adapters.md){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/y9c/cutseq){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Why CutSeq?

NGS library prep can need many sequential trimming and extraction steps. The
_SMARTer® Stranded Total RNA-Seq Kit v3_ workflow, for example, needs at least 9
operations — adapter removal, UMI extraction, linker masking, poly-A/T trimming —
done **in the correct order**. CutSeq compiles a single library *scheme* into one
native cutadapt pass, so the whole pipeline is one parameter and one command.

![](https://raw.githubusercontent.com/y9c/cutseq/main/docs/takaraV3.png)

## Key features

- **One-parameter trimming** — describe the whole library as a scheme and let
  CutSeq apply every step automatically.
- **Built-in + custom schemes** — 70+ bundled libraries (`cutseq --list-adapters`)
  or define your own grammar string.
- **Barcodes & UMIs** — capture and rename inline barcodes / UMIs into read names.
- **QC-aware** — discard reads with a reason (`too_short`, `too_many_n`, …).
- **Efficient & reproducible** — a single cutadapt pass, no intermediate I/O.

## Explore the docs

| Page | What it covers |
|------|----------------|
| [Quick Start](quickstart.md) | Install + first commands |
| [Adapter Schemes](adapters.md) | All built-in schemes, with an interactive viewer & `.gb`/`.dna` download |
| [Read Name Renaming](rename.md) | `--rename`, captures and transforms |
| [Spatial (DBiT-seq) libraries](spatial.md) | Parsing a spatial barcode-arm library |

---

{: .note }
CutSeq is open source and contributions are welcome — see the
[GitHub repo](https://github.com/y9c/cutseq).
