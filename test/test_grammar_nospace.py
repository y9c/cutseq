"""Unit tests for the space-free grammar tokenizer.

The scheme is written WITHOUT spaces. Tokens are split by character-class
transitions; adjacent same-class tokens (e.g. two lowercase inline barcodes)
merge into one. Numeric shorthand expands to the explicit run.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cutseq.grammar import parse_scheme, tokenize


def _describe(tokens):
    return [(t.kind, t.value) for t in tokens]


def test_adapter_then_capture_no_space():
    """Adapter run and N-capture abut without a space and split cleanly."""
    toks = tokenize("ACGTN8")
    assert _describe(toks) == [("adp", "ACGT"), ("capture", 8)]


def test_inline_adjacent_merge():
    """Adjacent lowercase inline barcodes merge into one token."""
    toks = tokenize("acgtcagt")
    assert _describe(toks) == [("inline", "acgtcagt")]


def test_mask_numeric_expands():
    """X3 collapses to a mask of length 3 (same as XXX)."""
    assert _describe(tokenize("X3")) == [("mask", 3)]
    assert _describe(tokenize("XXX")) == [("mask", 3)]


def test_capture_numeric_expands():
    """N5 collapses to a capture of length 5 (same as NNNNN)."""
    assert _describe(tokenize("N5")) == [("capture", 5)]
    assert _describe(tokenize("NNNNN")) == [("capture", 5)]


def test_full_scheme_no_space():
    """A full library layout written with no spaces tokenizes correctly."""
    # real: adapter + inline + capture, all abutted
    s = "ATCCACGTGCTTACGTcagtN8"
    toks = tokenize(s)
    assert _describe(toks) == [
        ("adp", "ATCCACGTGCTTACGT"),
        ("inline", "cagt"),
        ("capture", 8),
    ]


def test_insert_marker_orientation():
    """The + - : markers are extracted as single insert-marker tokens."""
    assert _describe(tokenize("ACGT+ACGT")) == [("adp", "ACGT"), ("insert", "+"), ("adp", "ACGT")]
    assert _describe(tokenize("ACGT-ACGT")) == [("adp", "ACGT"), ("insert", "-"), ("adp", "ACGT")]


def test_parse_scheme_no_space():
    """parse_scheme works on a space-free scheme, respecting the insert marker."""
    orientation, left, right = parse_scheme("ATCCACGTN8+AAGCAGTGG")
    assert orientation == "+"
    assert _describe(left) == [("adp", "ATCCACGT"), ("capture", 8)]
    assert _describe(right) == [("adp", "AAGCAGTGG")]


def test_polytail_after_adapter():
    """A poly tail directly after an adapter splits into adapter + polytail."""
    assert _describe(tokenize("ACGTAAA...A")) == [
        ("adp", "ACGT"),
        ("polytail", "A"),
    ]


def test_polytail_then_adapter():
    """A poly tail followed by another adapter parses cleanly."""
    assert _describe(tokenize("AAA...AAAGGG")) == [
        ("polytail", "A"),
        ("adp", "GGG"),
    ]


def test_polytail_does_not_swallow_insert():
    """An insert marker between an adapter and a poly tail is preserved."""
    assert _describe(tokenize("ACGT+TTT...TTT")) == [
        ("adp", "ACGT"),
        ("insert", "+"),
        ("polytail", "T"),
    ]


def test_homopolymer_adapter_numeric():
    """A12 expands to a 12-base homopolymer adapter."""
    assert _describe(tokenize("A12")) == [("adp", "AAAAAAAAAAAA")]
    assert _describe(tokenize("T10ACGT")) == [("adp", "TTTTTTTTTT"), ("adp", "ACGT")]


if __name__ == "__main__":
    import traceback

    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception:
                failures += 1
                print(f"FAIL {name}")
                traceback.print_exc()
    sys.exit(1 if failures else 0)
