"""Tests for inline-barcode auto-detection and the sequencing-primer database.

Covers:
  - terminal-fragment primer matching (partial, RC-aware)
  - splitting a barcode merged with p7 (3' side) / p5 (5' side)
  - no false positives on correct lowercase or known-primer-only schemes
  - custom schemes are never altered
  - --no-auto-inline and --list-primers CLI flags
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")
R1 = str(ROOT / "test" / "input_R1.fq.gz")
R2 = str(ROOT / "test" / "input_R2.fq.gz")

sys.path.insert(0, str(ROOT))
from cutseq import grammar  # noqa: E402
from cutseq.primers import is_known_primer, primer_name  # noqa: E402


def _kinds(scheme):
    _, left, right = grammar.parse_scheme(scheme)
    return [(t.kind, t.value) for t in left], [(t.kind, t.value) for t in right]


# --- primer matching -------------------------------------------------------


def test_known_primer_full_and_fragment():
    # full TruSeq R1 primer
    assert is_known_primer("ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
    # 3' read-adjacent fragment (as appears in schemes)
    assert is_known_primer("ACACGACGCTCTTCCGATCT")
    # truncated to the terminal-match minimum
    assert is_known_primer("GACGCTCTTCCGATCT")  # 15 nt
    # below the threshold is NOT a primer match
    assert not is_known_primer("GCTCTTCCGATCT")  # 12 nt


def test_known_primer_rc_and_mgi():
    # p7 in schemes is the RC of the R2 primer's 3' tail
    assert is_known_primer("AGATCGGAAGAGCACACGTC")
    # BGI / MGI
    assert is_known_primer("AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA")
    assert is_known_primer("AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG")


def test_random_sequence_not_a_primer():
    assert not is_known_primer("GGCCTAGGATCCGATCGTAG")


def test_primer_name_reports_match():
    assert "TruSeq p7 (read)" in primer_name("AGATCGGAAGAGCACACGTC")


# --- auto-detection --------------------------------------------------------


def test_split_barcode_merged_with_p7_right_side():
    # INLINE scheme where the user wrote ATCACG in uppercase next to p7
    left, right = _kinds(
        "AGTTCTACAGTCCGACGATCNNNNN+NNNNNATCACGAGATCGGAAGAGCACACGTC"
    )
    assert ("inline", "ATCACG") in right
    assert ("adp", "AGATCGGAAGAGCACACGTC") in right


def test_correct_lowercase_unchanged():
    left, right = _kinds(
        "AGTTCTACAGTCCGACGATCNNNNN+NNNNNatcacgAGATCGGAAGAGCACACGTC"
    )
    assert ("inline", "atcacg") in right
    assert ("adp", "AGATCGGAAGAGCACACGTC") in right


def test_split_barcode_merged_with_p5_left_side():
    left, right = _kinds("AGTTCTACAGTCCGACGATCGGATCC+AGATCGGAAGAGCACACGTC")
    assert ("inline", "GGATCC") in left
    assert ("adp", "AGTTCTACAGTCCGACGATC") in left


def test_custom_scheme_untouched():
    left, right = _kinds("GGGGCCCCAAAA+TTTTGGGGCCCCAAAA")
    assert left == [("adp", "GGGGCCCCAAAA")]
    assert right == [("adp", "TTTTGGGGCCCCAAAA")]


def test_builtin_scheme_untouched():
    # ECLIP10 outer adapters are known primers but there is no barcode to
    # detect, so output must be unchanged.
    left, right = _kinds("ACACGACGCTCTTCCGATCTXX-XNNNNNNNNNNAGATCGGAAGAGCACACGTC")
    assert ("adp", "ACACGACGCTCTTCCGATCT") in left
    assert ("adp", "AGATCGGAAGAGCACACGTC") in right


def test_auto_inline_can_be_disabled():
    _, _, right = grammar.parse_scheme(
        "AGTTCTACAGTCCGACGATCNNNNN+NNNNNATCACGAGATCGGAAGAGCACACGTC",
        auto_inline=False,
    )
    # no split happened; the merged run stays as one adapter
    assert any(t.kind == "adp" and t.value.startswith("ATCACG")
               for t in right)


def test_cli_no_auto_inline_flag():
    p = subprocess.run([CUTSEQ, "--no-auto-inline", "-A", "SMALLRNA",
                        "-n", R1, R2], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    # SMALLRNA has no barcode; flag should just be accepted and no-op
    assert "Auto-detected" not in p.stderr + p.stdout


def test_cli_list_primers():
    p = subprocess.run([CUTSEQ, "--list-primers"], capture_output=True, text=True)
    assert p.returncode == 0
    assert "TruSeq" in p.stdout and "Nextera" in p.stdout and "MGI" in p.stdout


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
