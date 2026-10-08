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


# --- genuine sequencing-site records ---------------------------------------
#
# The flat ``SEQUENCING_PRIMERS`` DB mixes flowcell anchors, ligation adapters
# and the *read-adjacent* sequencing primers. Only the last group are genuine
# sequencing-SITE records: they anneal upstream of where each read actually
# starts, so the read does not contain them. We record them here with an
# explicit role (R1 = read-1 site on the top strand, R2 = read-2 site on the
# top strand) and resolve them end-directed + orientation-aware, so a scheme's
# outer adapter can be recognised as a read boundary rather than in-read
# scaffold, while ligation / flowcell adapters stay read-visible.

# name -> (role, 5'->3' TOP-STRAND site).  A site occupies the OUTER position
# on its side of the scheme molecular map.  role 'R1' is the read-1 primer
# binding region at the left (5') edge; role 'R2' is the read-2 primer binding
# region at the right (3') edge.  The ACTUAL read-2 sequencing primer oligo is
# ``rc(site)`` (R2 reads the bottom strand), and an explicitly supplied
# ``--r2-primer`` is that oligo 5'->3'.
#
# Sources: support-docs.illumina.com/SHARE/AdapterSequences.  sRNA RA5/RA3 and
# MGI flowcell filters are NOT genuine read-adjacent sequencing SITES (they are
# ligation / flowcell adapters), so they are deliberately excluded here.
SEQUENCING_SITES = {
    # --- Illumina read-adjacent sequencing primers (genuine sites) ---
    # R1: the top-strand read-1 primer binding region.  A scheme's left outer
    # fragment must reach the primer's 3' EXTENSION END (a suffix of the full
    # oligo) to define a boundary.
    "TruSeq R1 (5')": ("R1", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"),
    "TruSeq p5 (read)": ("R1", "ACACGACGCTCTTCCGATCT"),
    # R2: the top-strand read-2 binding site.  A scheme's right outer fragment is
    # a PREFIX of the full top-strand site (the read-2 primer's 3' extension end
    # sits at the site's LEFT edge, i.e. the insert-adjacent end).
    "TruSeq R2 (3')": ("R2", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
    "TruSeq p7 (read)": ("R2", "AGATCGGAAGAGCACACGTC"),
    "Nextera R2": ("R2", "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"),
    "Nextera read (3')": ("R2", "CTGTCTCTTATACACATCT"),
}


def _site_boundary(candidate, site, side):
    """Return the read-start boundary offset within *candidate*, or ``None``.

    A genuine sequencing SITE defines where its read STARTS (the 3' extension
    end of the bound primer).  This is direction-specific:

    * R1 (left, top strand): the primer anneals at the read-1 binding region and
      extends toward the insert, so the read starts at the SITE's 3' end.  A
      scheme fragment is a boundary only if it reaches that 3' end (it is a
      SUFFIX of the full site, or the full site is a suffix of it, or the
      fragment equals a terminal piece ending at the site's 3' end).  The
      boundary offset is where the site's 3' end sits within *candidate*.
    * R2 (right, top strand): the read-2 primer's 3' extension end is at the
      SITE's LEFT (insert-adjacent) edge, so the read starts BEFORE the site.  A
      scheme fragment is a boundary only when it is a PREFIX of the full site;
      the boundary offset is at the start of the site within *candidate*.

    Same-boundary aliases may produce an equal offset; a candidate that matches
    only ``rc(site)`` or an interior piece is NOT a boundary.  Returns the
    within-candidate offset or ``None``.
    """
    c, p = _norm(candidate), _norm(site)
    if side == "R1":
        # The read starts at the site's 3' EXTENSION END.  A candidate is a
        # boundary iff its 3' end coincides with the site's 3' end: candidate is
        # the site, candidate ends with the site (site is a prefix), or candidate
        # is a 3'-terminal stem of the site (site ends with candidate).  In every
        # case the boundary offset is len(candidate) -- the read starts at the
        # candidate's 3' end.  (Interior fragments are rejected.)
        if c == p:
            return len(c)
        if len(c) > len(p) and c.endswith(p):
            return len(c)
        if len(c) < len(p) and p.endswith(c):
            return len(c)
        return None
    # R2: boundary at the site's left (insert-adjacent) edge.  A candidate is a
    # boundary iff its 5' end coincides with the site's 5' end: candidate is the
    # site, candidate is a PREFIX of the site, or candidate starts with the site.
    # The read starts at offset 0 of the site region.
    if c == p:
        return 0
    if p.startswith(c):
        return 0
    if len(c) > len(p) and c.startswith(p):
        return 0
    return None


def _site_boundaries(candidate, site, side):
    """Return ALL within-``candidate`` read-start boundary offsets for *site*.

    Unlike ``_site_boundary`` (which only accepts a terminal match), this also
    collects occurrences of the site EMBEDDED inside a merged token so a primer
    repeated at several offsets is detected (and can be reported as ambiguous)
    while a genuine single occurrence still resolves.  Returns a sorted de-duped
    list; empty when the candidate is not a boundary at all.

    * R1 (left, top strand): boundary at the site's 3' EXTENSION END.  Every full
      occurrence of the site at position ``p`` yields ``p + len(site)``; a
      candidate that is a 3'-terminal stem of the site (the site ends with it, or
      it equals a site whose 3' end reaches ``len(c)``) yields ``len(c)``.
    * R2 (right, top strand): boundary at the site's LEFT (insert-adjacent) edge.
      Every full occurrence at ``p`` yields ``p``; a candidate that is a prefix-
      stem of the site (the site starts with it, or it starts with the site)
      yields ``0``.
    """
    c, p = _norm(candidate), _norm(site)
    offs = []
    if side == "R1":
        start = 0
        while True:
            pos = c.find(p, start)
            if pos == -1:
                break
            offs.append(pos + len(p))
            start = pos + 1  # +1 catches overlapping (non-gapped) hits
        if c == p or (len(c) > len(p) and c.endswith(p)) or \
                (len(c) < len(p) and p.endswith(c)):
            offs.append(len(c))
    else:
        start = 0
        while True:
            pos = c.find(p, start)
            if pos == -1:
                break
            offs.append(pos)
            start = pos + 1
        if c == p or (len(c) > len(p) and c.startswith(p)) or \
                (len(c) < len(p) and p.startswith(c)):
            offs.append(0)
    return sorted(set(offs))


def sequencing_site_role(seq, side):
    """Return ``'R1'``/``'R2'`` if *seq* is a genuine sequencing-site record on
    the given ``side`` (``'R1'``/left or ``'R2'``/right top-strand position),
    else ``None``. Direction-aware and extension-end-anchored; ``rc(site)`` and
    interior pieces are rejected."""
    s = _norm(seq)
    for _name, (role, site) in SEQUENCING_SITES.items():
        if role != side:
            continue
        if _site_boundary(s, site, side) is not None:
            return role
    return None


def site_boundary_offset(seq, side):
    """The within-``seq`` read-start boundary offset for a recognised site on
    ``side``, else ``None``. This is what the read-plan uses to split a token:
    everything at/after the offset on the read's own arm is read-visible, the
    prefix (upstream) is skipped on THAT read."""
    s = _norm(seq)
    for _name, (role, site) in SEQUENCING_SITES.items():
        if role != side:
            continue
        off = _site_boundary(s, site, side)
        if off is not None:
            return off
    return None


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
