---
title: Home
layout: home
nav_order: 1
description: "CutSeq: Efficient, flexible adapter trimming for NGS data."
permalink: /
---

# ✂️ CutSeq
{: .fs-9 }

Trim sequencing adapters in a correct way.
{: .fs-6 .fw-300 }

[Get started](quickstart.md){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/y9c/cutseq){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Why CutSeq?

NGS library preparation can require many sequential trimming and extraction steps. For example, the _SMARTer® Stranded Total RNA-Seq Kit v3_ workflow needs at least 9 operations, including adapter removal, UMI extraction, linker masking, and poly-A/T trimming. CutSeq automates these steps in the correct order, reducing errors and saving time.

![](https://raw.githubusercontent.com/y9c/cutseq/main/docs/takaraV3.png)

## Key Features
- **Single-parameter trimming:** Specify your library scheme with one parameter.
- **Built-in and custom adapters:** Use a preset or define your own scheme.
- **Efficient and reproducible:** Avoids multiple I/O steps and manual errors.

---

{: .note }
CutSeq is open source and welcomes contributions! See the [GitHub repo](https://github.com/y9c/cutseq) for details.

## Quick Navigation

- [Quick Start](quickstart.md) - Installation and usage examples
- [Adapter Schemes](adapters.md) - Comprehensive guide to supported adapter patterns
