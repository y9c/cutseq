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

Or specify a custom adapter scheme:

```bash
cutseq -a "ACACGACGCTCTTCCGATCTXXX<XXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC" test_R1.fq.gz test_R2.fq.gz
```

![](https://raw.githubusercontent.com/y9c/cutseq/main/docs/explain_library.png)

The scheme string describes the structure of your library:
- Outermost: Illumina adapters
- Inner: Inline barcodes or custom primers (in parentheses)
- Next: UMI sequences (`N`)
- Next: Masked sequences (`X`)
- Center: Library insert (`>`, `<`, or `-` for orientation)

See the [documentation](https://cutseq.yech.science) for more details.
