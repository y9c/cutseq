"""Regression tests for grammar-engine edge-case bugs:
  - zero-length capture/mask/adapter numeric shorthand (N0 / X0 / A0) must be
    rejected instead of silently truncating reads to empty
  - a defensive guard in _ConditionalCutter (length == 0) never wipes a read
  - a '>' back adapter on the R2 side is emitted in single-end mode too
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dnaio import SequenceRecord  # noqa: E402
from cutadapt.info import ModificationInfo  # noqa: E402

from cutseq.grammar import (  # noqa: E402
    CompiledScheme,
    _ConditionalCutter,
    parse_scheme,
)


def _read(seq, name="r"):
    return SequenceRecord(name, seq, "I" * len(seq))


def test_zero_length_numeric_shorthand_rejected():
    for bad in ("ACGT+N0", "ACGT+X0"):
        try:
            parse_scheme(bad)
        except ValueError as e:
            assert ">= 1" in str(e), e
        else:
            raise AssertionError(f"expected ValueError for {bad!r}")


def test_zero_length_adapter_shorthand_rejected():
    try:
        parse_scheme("A0+ACGT")
    except ValueError as e:
        assert ">= 1" in str(e), e
    else:
        raise AssertionError("expected ValueError for 'A0'")


def test_conditional_cutter_zero_length_is_noop():
    # Defensive guard: even a hand-built length-0 cutter must not wipe reads.
    read = _read("ACGTACGTACGT")
    mod = _ConditionalCutter(0)
    out = mod(read, ModificationInfo(read))
    assert out is read


def test_back_adapter_on_r2_side_emitted_in_single_end():
    # '>SEQ' on the right side (R2) must still trim the 3' end in single-end.
    o, l, r = parse_scheme("AGTTCTACAGTCCGACGATC+>AGATCGGAAGAGCACACGTC")
    cs = CompiledScheme(o, l, r)
    mods, r2 = cs.modifiers(paired=False)
    assert r2 == []
    # front adapter + the R2 back adapter mirrored to R1's 3' end
    assert len(mods) == 2
    assert type(mods[1]).__name__ == "AdapterCutter"
    assert "AGATCGGAAGAGCACACGTC" in mods[1].adapters[0].sequence


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
