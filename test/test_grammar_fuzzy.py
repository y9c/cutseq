"""Unit tests for the grammar engine's fuzzy adapter / inline matching.

The grammar engine uses cutadapt's default FUZZY matchers (max_errors=0.1,
min_overlap=3), tolerating mismatches within that error rate — not strict
exact prefix equality.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dnaio import SequenceRecord
from cutadapt.info import ModificationInfo

import cutadapt.modifiers as mods
from cutadapt.adapters import FrontAdapter, PrefixAdapter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cutseq"))

from cutseq.grammar import _ConditionalCutter  # noqa: E402


def _read(seq, name="r"):
    return SequenceRecord(name, seq, None)


def test_adapter_fuzzy_match_tolerates_one_mismatch():
    """A 5' adapter with a single mismatch must still be trimmed (fuzzy)."""
    # 11nt adapter: 1 mismatch / 11 = 0.09 < cutadapt default error rate 0.1
    mod = mods.AdapterCutter([FrontAdapter("AAGCAGTGGTA")])
    read = _read("AAGCAGTGGCA" + "AAAACCGGTT" * 3)  # 1 mismatch at pos 9
    info = ModificationInfo(read)
    out = mod(read, info)
    # adapter fully trimmed off -> sequence is the insert
    assert str(out.sequence) == "AAAACCGGTT" * 3
    assert len(info.matches) == 1


def test_inline_fuzzy_match_captures_matched_sequence():
    """An inline barcode with mismatches is matched and trimmed (fuzzy)."""
    # 11nt barcode: 1 mismatch is within the default 0.1 error rate
    mod = mods.AdapterCutter([PrefixAdapter("ACGTACGTACG")])
    read = _read("ACGTACGTATG" + "GGGGCCCC" * 3)  # 1 mismatch at pos 9
    info = ModificationInfo(read)
    out = mod(read, info)
    assert len(info.matches) == 1, "capture must be recorded"
    # matched region length must equal the barcode
    assert str(out.sequence) == "GGGGCCCC" * 3


def test_inline_no_match_untouched():
    """A non-matching inline barcode leaves the read untouched (no capture)."""
    mod = mods.AdapterCutter([PrefixAdapter("ACGTACGT")])
    read = _read("TTTTTTTT" + "GGGGCCCC" * 3)
    info = ModificationInfo(read)
    out = mod(read, info)
    assert str(out.sequence) == "TTTTTTTT" + "GGGGCCCC" * 3
    assert info.matches == []


def test_mask_still_exact():
    """Mask / fixed-length trims remain exact-length (no fuzzy change)."""
    mod = _ConditionalCutter(4, conditional=False)
    read = _read("NNNN" + "ACGTACGT")
    info = ModificationInfo(read)
    out = mod(read, info)
    assert str(out.sequence) == "ACGTACGT"


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
