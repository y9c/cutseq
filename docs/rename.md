---
title: Read Name Renaming
nav_order: 3
---

# Read Name Renaming

By default, CutSeq appends every captured UMI / inline barcode to the read
name with `_`, exactly like the legacy engine (`@READID_CTATTAAAAA`). The
`--rename` option (alias `--name-format`) replaces that with a custom
template, and — crucially — lets you place each **individual capture**
anywhere you want, with optional transforms.

```bash
cutseq -A INLINE --rename '{id}_BC1:{1}_BC2:{2}_umi:rc({3})' in_R1.fq.gz in_R2.fq.gz
```

{: .note }
The template string is parsed **once** when the pipeline is built and the
result is compiled into a lazy rename function. Each read then only walks the
compiled field list — there is no per-read string parsing.

## What the template can contain

The template is a mix of literal text, **variables** in `{...}`, and
**function calls**.

### Native cutadapt variables

`{id}`, `{header}`, `{comment}`, `{cut_prefix}`, `{cut_suffix}`,
`{adapter_name}`, `{match_sequence}`, `{rc}` — byte-identical to
`cutadapt --rename`.

### Positional captures

`{1}`, `{2}`, ... reference the captured parts in **scheme written order**
(left side first, then right side), 1-based. Both `N`-captures (UMIs) and
inline barcodes are counted. For example, in
`AGTTCTACAGTCCGACGATCNNNNN+NNNNNatcacgAGATCGGAAGAGCACACGTC`:

| capture | part | anchor read |
|---|---|---|
| `{1}` | left `NNNNN` UMI | R1 |
| `{2}` | right `NNNNN` UMI | R2 (read as its reverse complement) |
| `{3}` | inline `atcacg` | R2 |

Captures can also be referenced by a `label` (from a YAML scheme `label:`).

### Transform functions

Any capture may be wrapped in functions. They **nest** and are
**case-insensitive** (both `rc` and `RC` work):

| function | meaning | example → result (`x = "AACC"`) |
|---|---|---|
| `rc(x)` | reverse complement | `rc(x)` → `GGTT` |
| `rev(x)` | reverse only | `rev(x)` → `CCAA` |
| `comp(x)` | complement only | `comp(x)` → `TTGG` |
| `canon(x)` | canonical UMI: `min(x, rc(x))` | `canon(x)` → `AACC` |
| `upper(x)` | uppercase | `upper({1})` |
| `lower(x)` | lowercase | `lower({1})` |
| `len(x)` | length | `len(x)` → `4` |
| `left(x, k)` | first `k` bases (`k<0`: drop last \|k\|) | `left(x,2)` → `AA` |
| `right(x, k)` | last `k` bases (`k<0`: drop first \|k\|) | `right(x,2)` → `CC` |
| `slice(x, a, b)` | 1-based **inclusive** range (negatives allowed) | `slice(x,1,2)` → `AA` |

Nesting example: `upper(RC({1}))`, `rc(comp({2}))`, `slice({1},1,4)`.

## Paired-end semantics

- An **unprefixed** capture resolves to its *anchor* read
  (left-side captures → R1, right-side → R2; see the table above), so both
  mates get the **same** extracted value. Left-side captures are read in R1
  orientation; right-side captures are read in R2 orientation (i.e. the
  reverse complement of the top-strand sequence) — matching legacy naming.
- `{r1.1}` / `{r2.1}` force a specific mate; each `info` also retains the
  *mirrored* value of the other side's capture (e.g. `{r2.1}` is the reverse
  complement of `{r1.1}`).

```bash
# same value on both mates:
cutseq -A INLINE --rename '{id}_BC1:{1}_BC2:{2}' in_R1.fq.gz in_R2.fq.gz
# force sides explicitly:
cutseq -A INLINE --rename '{id}_fwd={r1.1}_rev={r2.1}' in_R1.fq.gz in_R2.fq.gz
```

## Examples

```bash
# Label barcodes and reverse-complement the inline barcode (INLINE: 3 parts)
cutseq -A INLINE --rename '{id}_BC1:{1}_BC2:{2}_umi:rc({3})' \
    in_R1.fq.gz in_R2.fq.gz
# @READID_BC1:GGGNC_BC2:CNCCC_umi:ATCACG

# Canonical UMI (min(umi, rc(umi))) for UMI collapsing (ECLIP10: 1 capture)
cutseq -A ECLIP10 --rename '{id}_c=canon({1})' in_R1.fq.gz

# Truncate a long UMI to 6 bases and uppercase the reverse complement
# (TAKARAV3: 8-nt UMI, 1 capture)
cutseq -A TAKARAV3 --rename '{id}_{left(rc({1}),6)}' in_R1.fq.gz in_R2.fq.gz

# All captures appended, with ':' between them
# (reproduces the legacy --capture-separator ':' behavior in --rename)
cutseq -A INLINE --rename '{id}:{1}:{2}' in_R1.fq.gz in_R2.fq.gz
# @READID:GGGNC:CNCCC
```

{: .note }
Templates that contain **no** capture references are handed to cutadapt's own
renamer unchanged (byte-identical behavior). The internal engine only runs
when `{1}`… `{N}`, labels, or functions are used, so the fast default path is
never slowed down.

## Errors

Unknown variables, out-of-range captures (`{4}` when the scheme has three
parts), missing numeric arguments (`left({1})`), and unbalanced parentheses
all produce a clear error message at build time — before any reads are
processed.

## See also

- [Quick start](quickstart.md)
- [Adapter schemes](adapters.md)
