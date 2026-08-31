"""Known Illumina / BGI (MGI) sequencing primer and adapter sequences.

Used to auto-detect inline barcodes in a library scheme. A scheme is expected
to carry the sequencing primers (p5/p7 or the read-adjacent adapter sequences)
at its two outermost ends; any fixed uppercase sequence between them is an
inline barcode (matched and trimmed), so a barcode written in uppercase by
[Output truncated. Continue viewing below]
mistake can be detected and treated as ``inline``.

Sources:
  - Illumina Adapter Sequences document 1000000002694 (support.illumina.com)
  - teichlab.github.io/scg_lib_structs/methods_html/Illumina.html
  - MGI / DNBSEQ oligo documentation (MGIEasy UDB, NEBNext Multiplex for MGI)
  - OpenGene/fastp issue #259 (MGI/BGI adapter sequences)

Sequences are stored 5' -> 3' in the top-strand orientation.
"""

import logging

# name -> 5'->3' sequence. Grouped by platform / library type.
SEQUENCING_PRIMERS = {
    # --- Illumina flowcell anchors (full oligos) ---
    "P5 (flowcell)": "AATGATACGGCGACCACCGAGATCTACAC",
    "P7 (flowcell)": "CAAGCAGAAGACGGCATACGAGAT",

    # --- TruSeq read-adjacent adapters (appear at read ends) ---
    "TruSeq R1 (5')": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "TruSeq R2 (3')": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "TruSeq p5 (read)": "ACACGACGCTCTTCCGATCT",
    "TruSeq p7 (read)": "AGATCGGAAGAGCACACGTC",
    "TruSeq universal adapter": "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "TruSeq 3' adapter": "AGATCGGAAGAGCACACGTCT",
    "TruSeq RiboProfile fwd primer": "ATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACG",

    # --- TruSeq small RNA ---
    "TruSeq sRNA RA5 (5')": "GTTCAGAGTTCTACAGTCCGACGATC",
    "TruSeq sRNA RA3 (3')": "TGGAATTCTCGGGTGCCAAGG",
    "TruSeq sRNA RT primer": "GCCTTGGCACCCGAGAATTCCA",
    "TruSeq sRNA RP1": "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",

    # --- Nextera transposase / read adapters ---
    "Nextera R1": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "Nextera R2": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "Nextera read (5')": "AGATGTGTATAAGAGACAG",
    "Nextera read (3')": "CTGTCTCTTATACACATCT",

    # --- BGI / MGI / DNBSEQ ---
    "MGI fwd": "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    "MGI rev": "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
    "MGI universal": "AAGTCGGA",
}


def _norm(seq):
    return str(seq).upper().replace("U", "T").translate(
        str.maketrans("", "", " *:+-.'\u2032\u00b0")
    )


def _rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


_KNOWN = {_norm(s) for s in SEQUENCING_PRIMERS.values()}

# A primer fragment in a scheme is the insert-adjacent (terminal) portion of
# the full oligo, so we accept terminal matches down to this length. A 15-mer
# is ~4^15 ≈ 1e9 combinations — effectively unique against this small,
# curated primer database.
MIN_PRIMER_MATCH = 15


def _terminal_match(candidate, primer):
    """True if *candidate* and *primer* share a terminal fragment of at least
    ``MIN_PRIMER_MATCH`` bp. The scheme adapter may be either shorter than the
    full oligo (the common case) or, for very short oligos, longer than it."""
    c, p = _norm(candidate), _norm(primer)
    if len(c) >= MIN_PRIMER_MATCH and (p.endswith(c) or p.startswith(c)):
        return True
    if len(p) >= MIN_PRIMER_MATCH and (c.endswith(p) or c.startswith(p)):
        return True
    return False


def _matches_any(candidate):
    for primer in _KNOWN:
        if _terminal_match(candidate, primer) or _terminal_match(candidate, _rc(primer)):
            return primer
    return None


def is_known_primer(seq):
    """True if *seq* matches a known sequencing primer (terminal fragment,
    either strand, down to ``MIN_PRIMER_MATCH`` bp)."""
    return _matches_any(_norm(seq)) is not None


def primer_name(seq):
    """Return the name(s) of a matched primer, else None."""
    s = _norm(seq)
    names = []
    for name, p in SEQUENCING_PRIMERS.items():
        pn = _norm(p)
        if _terminal_match(s, pn) or _terminal_match(s, _rc(pn)):
            names.append(name)
    return names or None


def detect_5prime_primer(seq, max_scan=40):
    """Check whether a read's 5' end begins with a known sequencing primer.

    The read may start at the primer's 5' end (full oligo read-through) or
    deeper in (the 5' part was clipped); any terminal fragment of at least
    ``MIN_PRIMER_MATCH`` bp on either strand counts. Returns
    ``(fragment_len, primer_name, matched_fragment)`` for the longest match,
    or ``None``.
    """
    s = _norm(seq)[:max_scan]
    if len(s) < MIN_PRIMER_MATCH:
        return None
    best = None  # (fragment_len, name, fragment)
    for name, p in SEQUENCING_PRIMERS.items():
        pn = _norm(p)
        for frag in (pn, _rc(pn)):
            if len(frag) < MIN_PRIMER_MATCH:
                continue
            limit = min(len(s), len(frag))
            L = limit
            while L >= MIN_PRIMER_MATCH:
                if s[:L] == frag[:L]:
                    if best is None or L > best[0]:
                        best = (L, name, frag)
                    break
                L -= 1
    return best


def detect_5prime_from_reads(paths, n=200, max_scan=40):
    """Detect the read-5' sequencing primer(s) from a sample of reads.

    ``paths`` is one (R1) or two (R1, R2) input FASTQ paths. Returns a list
    aligned with the inputs: for each file, ``(primer_name_or_None,
    primer_seq_or_None)`` where ``primer_seq`` is the matched fragment to trim
    (most common detection across the sample). Ungzipped/gzipped reads are both
    supported.
    """
    try:
        import dnaio
    except ImportError:  # pragma: no cover
        logging.debug("dnaio unavailable; primer auto-detection skipped")
        return [(None, None) for _ in paths]

    out = []
    for path in paths:
        counts = {}      # (name, frag) -> count
        read_total = 0
        with dnaio.open(path, fileformat="fastq") as fh:
            for rec in fh:
                read_total += 1
                hit = detect_5prime_primer(rec.sequence, max_scan)
                if hit is not None:
                    _L, name, frag = hit
                    counts[(name, frag)] = counts.get((name, frag), 0) + 1
                if read_total >= n:
                    break
        if not counts:
            out.append((None, None))
            continue
        (name, frag), _ = max(counts.items(), key=lambda kv: kv[1])
        out.append((name, frag))
    return out
