"""Focused tests for the DEFAULT (no ``--rename``) capture-naming contract.

These cover the behaviour the read-plan owner left to the naming lane: the
default concatenates EVERY captured segment (``N`` UMIs + inline barcodes) in
SCHEME WRITTEN ORDER — not the R2-execution ``cut_prefix``/``cut_suffix`` order
the legacy fast path used — and obeys an all-or-nothing completeness policy
(any required anchor capture missing/incomplete -> the WHOLE combined suffix is
empty, ``id_``), never a partial plausible concat. It also pins the duplicate
label rules (across BOTH arms and auto-generated collisions) and the per-read
registry cleanliness across sequential renamers.

Completeness is judged by presence of the anchored capture in the per-read
registry (``_record_capture`` is only called for the segments the walk actually
read out; a capture the anchor read did not reach is simply absent).

These tests drive the same public pipeline the CLI uses (``_build_scheme`` +
``_scheme_modifiers``) with the DEFAULT renamer (``name_format=None``).
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cutseq import grammar
from cutseq.common import BUILDIN_ADAPTERS
from cutseq.run import CutadaptConfig, _build_scheme, _scheme_modifiers
from dnaio import SequenceRecord
from cutadapt.info import ModificationInfo


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# --- INLINE scheme + controlled molecule (same molecule as rename-engine) ----

_SCHEME = BUILDIN_ADAPTERS["INLINE"]
_P5, _P7 = "AGTTCTACAGTCCGACGATC", "AGATCGGAAGAGCACACGTC"
_INS = "T" * 40
# Top strand 5' -> 3': p5 | umi5 | insert | umi3 | bc(ATCACG) | p7
_R1_FULL = _P5 + "AAAAA" + _INS + "CCCCC" + "ATCACG" + _P7
_R2_FULL = _rc("ATCACG") + _rc("CCCCC") + _rc(_INS) + _rc("AAAAA") + _rc(_P5)
# Written capture order (values as read on the ANCHOR mate):
#   capture1 = umi5  (anchor R1) -> AAAAA
#   capture2 = umi3  (anchor R2) -> rc(CCCCC)      = GGGGG
#   capture3 = inline (anchor R2) -> rc(ATCACG)    = CGTGAT
_PE_DEFAULT = "AAAAA" + "GGGGG" + "CGTGAT"


# --- single-end sample scheme -----------------------------------------------

_SE_SCHEME = "ACGTACGTN4GTACN4"
_SE_FULL = "ACGTACGT" + "AACC" + "GTAC" + "GGTT" + "CCCGGGAtc"


# --- helpers -----------------------------------------------------------------


def _run_default_pair(r1_seq, r2_seq, scheme=_SCHEME):
    """Run the compiled PAIRED modifiers with the DEFAULT renamer (no --rename)."""
    s = CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=True, settings=s)  # mods[-1] = default renamer
    r1 = SequenceRecord("x/1", r1_seq, "I" * len(r1_seq))
    r2 = SequenceRecord("x/2", r2_seq, "I" * len(r2_seq))
    i1, i2 = ModificationInfo(r1), ModificationInfo(r2)
    try:
        for step in mods:
            if isinstance(step, tuple):
                m1, m2 = step
                n1 = m1(r1, i1) if m1 else r1
                n2 = m2(r2, i2) if m2 else r2
                r1, r2 = n1 or r1, n2 or r2
            else:
                res = step(r1, r2, i1, i2)
                if res is not None:
                    r1, r2 = res
        return r1, r2
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def _run_default_single(read_seq, scheme=_SE_SCHEME):
    """Run the compiled SINGLE-end modifiers with the DEFAULT renamer."""
    s = CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    read = SequenceRecord("x", read_seq, "I" * len(read_seq))
    info = ModificationInfo(read)
    try:
        for m in mods:
            out = m(read, info)
            if out is not None:
                read = out
        return read
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def _default_renamer(scheme, paired):
    """Build the default renamer for a scheme string (via the public path)."""
    orientation, left, right = grammar.parse_scheme(scheme)
    cs = grammar.CompiledScheme(orientation, left, right)
    return cs.renamer(paired=paired), left, right


# --- default: written-order concatenation (no --rename) ----------------------


def test_default_paired_written_order_all_segments():
    """The INLINE scheme has captures on BOTH arms; the DEFAULT name must join
    umi5 + umi3 + inline in SCHEME WRITTEN order on both mates."""
    r1, r2 = _run_default_pair(_R1_FULL, _R2_FULL)
    assert r1.sequence == _INS, r1.sequence
    assert r2.sequence == _rc(_INS), r2.sequence
    assert r1.name == f"x/1_{_PE_DEFAULT}", r1.name
    assert r2.name == f"x/2_{_PE_DEFAULT}", r2.name


def test_default_paired_no_mirror_duplicate():
    """In a FULL read-through both mates observe every capture (R1 its 3' mirror
    of umi3/inline, R2 its 5' own arm). The default must take each logical
    capture ONCE from its anchor read, never double-concatenate the mirrored
    values."""
    r1, r2 = _run_default_pair(_R1_FULL, _R2_FULL)
    # If the renamer naively merged both mates' registries the name would be
    # doubled (A... + G... + C... twice). Exact equality proves it is once.
    assert r1.name.count("AAAAA") == 1
    assert r1.name.count("GGGGG") == 1
    assert r1.name.count("CGTGAT") == 1
    assert r1.name == f"x/1_{_PE_DEFAULT}"
    assert r2.name == f"x/2_{_PE_DEFAULT}"


def test_default_single_written_order_all_segments():
    """Single-end, two adjacent N captures: default joins them in written order
    (first then second), not read-processing order."""
    read = _run_default_single(_SE_FULL)
    assert read.sequence == "CCCGGGAtc", read.sequence
    assert read.name == "x_AACCGGTT", read.name


def test_default_scheme_without_captures_keeps_plain_id():
    """A scheme with no capture/inline parts keeps the plain ``{id}`` (no
    trailing separator) under the default renamer."""
    read = _run_default_single("ACGTACGT" + "GGGGGGGG", scheme="ACGTACGT")
    assert read.name == "x", read.name


# --- default: all-or-nothing completeness -----------------------------------


def test_default_incomplete_empties_whole_composite():
    """Multisegment default: if ANY required capture is missing/incomplete the
    ENTIRE combined suffix is empty (``id_``) — never a partial plausible concat.
    Here the second of two SE captures is absent from the registry."""
    renamer, left, right = _default_renamer(_SE_SCHEME, paired=False)
    labels, _anchors = grammar._capture_meta(left, right)
    assert len(labels) == 2, labels
    read = SequenceRecord("x", "A" * 40, "I" * 40)
    info = ModificationInfo(read)
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        # only the FIRST capture is complete; the second is a required anchor
        # capture that was not read out (absent) -> whole suffix empty
        grammar._record_capture(info, labels[0], "AACC")
        renamer(read, info)
        assert read.name == "x_", read.name
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def test_default_multicapture_missing_single_anchor_also_empty():
    """Paired: with BOTH arms configured but only the RIGHT arm's captures
    complete (the LEFT/anchor-R1 capture absent), the whole composite empties."""
    orientation, left, right = grammar.parse_scheme(_SCHEME)
    cs = grammar.CompiledScheme(orientation, left, right)
    renamer = cs.renamer(paired=True)
    labels, anchors = grammar._capture_meta(left, right)
    # left (R1) barcode1 = umi5; right (R2) barcode2 = umi3, barcode3 = inline
    read1 = SequenceRecord("x/1", "A" * 20, "I" * 20)
    read2 = SequenceRecord("x/2", "T" * 20, "I" * 20)
    i1, i2 = ModificationInfo(read1), ModificationInfo(read2)
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        # complete ONLY the RIGHT-arm captures; the LEFT (anchor R1) is missing
        for label, anchor in zip(labels, anchors):
            if anchor == 1:
                grammar._record_capture(i2, label, "GGGGG" if label.endswith("2") else "CGTGAT")
        renamer(read1, read2, i1, i2)
        assert read1.name == "x/1_", read1.name
        assert read2.name == "x/2_", read2.name
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def test_default_complete_both_arms_still_joins():
    """Control: with EVERY anchor capture present the same PE scheme joins in
    written order (not left-empty because a single arm is populated)."""
    orientation, left, right = grammar.parse_scheme(_SCHEME)
    cs = grammar.CompiledScheme(orientation, left, right)
    renamer = cs.renamer(paired=True)
    labels, anchors = grammar._capture_meta(left, right)
    read1 = SequenceRecord("x/1", "A" * 20, "I" * 20)
    read2 = SequenceRecord("x/2", "T" * 20, "I" * 20)
    i1, i2 = ModificationInfo(read1), ModificationInfo(read2)
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        for label, anchor in zip(labels, anchors):
            src = i2 if anchor == 1 else i1
            grammar._record_capture(src, label, "AAAAA" if anchor == 0 else ("GGGGG" if label.endswith("2") else "CGTGAT"))
        renamer(read1, read2, i1, i2)
        assert read1.name == "x/1_AAAAAGGGGGCGTGAT", read1.name
        assert read2.name == "x/2_AAAAAGGGGGCGTGAT", read2.name
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


# --- duplicate labels (rejected, never silently overwritten) ----------------


def _write_yaml_parts(parts):
    import tempfile
    import yaml

    td = tempfile.TemporaryDirectory()
    path = f"{td.name}/s.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump({"parts": parts}, fh)
    return td, path


def test_duplicate_label_across_both_arms_rejected():
    """The SAME explicit label on a LEFT (R1) capture and a RIGHT (R2) capture is
    ambiguous and must fail actionably at scheme construction — not silently
    overwrite one segment with the other."""
    parts = [
        {"type": "capture", "length": 4, "label": "dup"},
        {"type": "insert", "value": "-"},
        {"type": "capture", "length": 4, "label": "dup"},
    ]
    td, path = _write_yaml_parts(parts)
    try:
        s = CutadaptConfig()
        with pytest.raises(ValueError):
            _build_scheme(path, s)
    finally:
        td.cleanup()


def test_duplicate_label_same_arm_rejected():
    """Two captures on the SAME arm with the same label fail actionably."""
    parts = [
        {"type": "capture", "length": 4, "label": "dup"},
        {"type": "adapter", "seq": "ACGT"},
        {"type": "capture", "length": 4, "label": "dup"},
    ]
    td, path = _write_yaml_parts(parts)
    try:
        s = CutadaptConfig()
        with pytest.raises(ValueError):
            _build_scheme(path, s)
    finally:
        td.cleanup()


def test_auto_barcode_collision_with_explicit_label_rejected():
    """An UNLABELED capture auto-targeted at ``barcode{n}`` that collides with an
    EXPLICIT ``barcode{n}`` label is a collision and must fail actionably."""
    parts = [
        {"type": "capture", "length": 4},                    # auto -> barcode1
        {"type": "insert", "value": "-"},
        {"type": "capture", "length": 4, "label": "barcode1"},  # collision
    ]
    td, path = _write_yaml_parts(parts)
    try:
        s = CutadaptConfig()
        with pytest.raises(ValueError):
            _build_scheme(path, s)
    finally:
        td.cleanup()


# --- registry cleanliness across sequential renamers ------------------------


def test_sequential_default_renamers_no_registry_leak():
    """Each default renamer consumes and pops the per-read registry entry, so
    building/running renamers back-to-back must not leak entries between reads
    or leave the module-global registry dirty."""
    grammar._RENAME_NEEDS_CAPTURES = False
    grammar._capture_registry.clear()
    for scheme in (_SE_SCHEME, _SE_SCHEME):
        orientation, left, right = grammar.parse_scheme(scheme)
        cs = grammar.CompiledScheme(orientation, left, right)
        ren = cs.renamer(paired=False)
        labels, _anchors = grammar._capture_meta(left, right)
        read = SequenceRecord("x", "A" * 30, "I" * 30)
        info = ModificationInfo(read)
        for label in labels:
            grammar._record_capture(info, label, "AAAA")
        ren(read, info)
        # consumed for THIS read -> nothing left behind for the next read
        assert grammar._capture_registry == {}, grammar._capture_registry
    grammar._RENAME_NEEDS_CAPTURES = False
    grammar._capture_registry.clear()


def test_default_single_umi_preserved_single_capture():
    """A scheme with a single capture keeps the ``{id}_{value}`` legacy shape
    (not bumped to an empty composite)."""
    read = _run_default_single("ACGTACGT" + "AACC" + "GG", scheme="ACGTACGTN4")
    assert read.sequence == "GG", read.sequence
    assert read.name == "x_AACC", read.name