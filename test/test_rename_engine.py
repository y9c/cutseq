"""Tests for the labeled-capture --rename engine.

Uses the built-in INLINE scheme
``AGTTCTACAGTCCGACGATCNNNNN+NNNNNatcacgAGATCGGAAGAGCACACGTC`` and a controlled
synthetic molecule::

    p5 | umi5(AAAAA) | insert | umi3(CCCCC) | bc(ATCACG) | p7

Capture order (written): capture1 = umi5 (anchor R1, value AAAAA),
capture2 = umi3  (anchor R2, value rc(CCCCC)=GGGGG, read-oriented),
capture3 = inline atcacg (anchor R2, value rc(ATCACG)=CGTGAT).

Covers positional captures, transforms (rc/rev/comp/canon/len/left/right/
slice/upper/lower), nesting, paired side-forcing (r1./r2.), identical names
across mates, and clean errors.
"""

import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")

sys.path.insert(0, str(ROOT))
from cutseq import grammar  # noqa: E402
from cutseq.common import BUILDIN_ADAPTERS  # noqa: E402
from cutseq.run import CutadaptConfig, _build_scheme, _scheme_modifiers  # noqa: E402

from dnaio import SequenceRecord  # noqa: E402
from cutadapt.info import ModificationInfo  # noqa: E402


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


_SCHEME = BUILDIN_ADAPTERS["INLINE"]
_P5, _P7 = "AGTTCTACAGTCCGACGATC", "AGATCGGAAGAGCACACGTC"
_INS = "T" * 40
_R1 = _P5 + "AAAAA" + _INS + "CCCCC" + "ATCACG" + _P7
_R2 = _rc(_P7) + _rc("ATCACG") + _rc("CCCCC") + _rc(_INS) + _rc("AAAAA") + _rc(_P5)


def _run_paired(name_format):
    """Run the compiled paired modifiers on one synthetic pair; return names."""
    s = CutadaptConfig()
    cs = _build_scheme(_SCHEME, s)
    mods = _scheme_modifiers(cs, paired=True, settings=s)
    mods[-1] = cs.renamer(paired=True, name_format=name_format)
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        r1 = SequenceRecord("x/1", _R1, "I" * len(_R1))
        r2 = SequenceRecord("x/2", _R2, "I" * len(_R2))
        i1, i2 = ModificationInfo(r1), ModificationInfo(r2)
        for step in mods:
            if isinstance(step, tuple):
                m1, m2 = step
                n1 = m1(r1, i1) if m1 else r1
                n2 = m2(r2, i2) if m2 else r2
                r1, r2 = n1 or r1, n2 or r2
            else:
                step(r1, r2, i1, i2)
        return r1.name, r2.name
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def _name(format_):
    n1, _ = _run_paired(format_)
    return n1


# --- positional + transforms -------------------------------------------------


def test_positional_captures():
    n1 = _name("{id}_c1={1}_c2={2}_c3={3}")
    assert n1 == "x/1_c1=AAAAA_c2=GGGGG_c3=CGTGAT"
    # both mates get identical values for unprefixed captures
    assert _run_paired("c={1}") == ("c=AAAAA", "c=AAAAA")


def test_rc_rev_comp_canon():
    n1 = _name("{id}_rc=rc({1})_rev=rev({1})_comp=comp({1})_canon=canon({2})")
    # rc(AAAAA)=TTTTT ; rev(AAAAA)=AAAAA ; comp(AAAAA)=TTTTT
    assert "_rc=TTTTT_rev=AAAAA_comp=TTTTT" in n1
    # canon(capture2)=min(GGGGG, rc(GGGGG)=CCCCC)=CCCCC
    assert n1.endswith("_canon=CCCCC")


def test_upper_lower_len_left_right_slice():
    n1 = _name("{id}_u=upper(rc({1}))_l=len({2})_lf=left({2},2)"
               "_rt=right({2},2)_sl=slice({1},2,4)_lo=lower({1})")
    assert "_u=TTTTT_l=5_lf=GG_rt=GG_sl=AAA_lo=aaaaa" in n1


def test_nested_functions_and_case_insensitive_names():
    n1 = _name("{id}_n=upper(RC({1}))_d=rc(comp({2}))")
    assert "_n=TTTTT_d=GGGGG" in n1


def test_paired_side_forcing():
    n1, n2 = _run_paired("{id}_a={r1.1}_b={r2.1}_c={r2.2}_d={r1.2}")
    # r2.1 = R2 mirror of capture1 = rc(AAAAA) = TTTTT; r2.2 = GGGGG
    assert n1 == "x/1_a=AAAAA_b=TTTTT_c=GGGGG_d=CCCCC"
    assert n2 == "x/2_a=AAAAA_b=TTTTT_c=GGGGG_d=CCCCC"


# --- errors ------------------------------------------------------------------


def test_out_of_range_capture_errors():
    try:
        _name("{id}_{4}")
    except ValueError as e:
        assert "only" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_missing_numeric_argument_errors():
    try:
        _name("{id}_left({1})")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_unknown_variable_errors():
    try:
        # the transformed '.rc' suffix syntax is not supported; use rc({n})
        _name("{id}_{1.rc}")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


# --- CLI end-to-end ----------------------------------------------------------


def _cli(template, *extra):
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            f1.write(f"@x/1\n{_R1}\n+\n{'I' * len(_R1)}\n")
            f2.write(f"@x/2\n{_R2}\n+\n{'I' * len(_R2)}\n")
        p = subprocess.run(
            [CUTSEQ, "-A", "INLINE", "-O", f"{td}/o", "--rename", template,
             "-m", "10", *extra, f"{td}/t1.fq", f"{td}/t2.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        with gzip.open(f"{td}/o_trimmed_R1.fastq.gz", "rt") as f:
            return next(f).strip()


def test_cli_rename_end_to_end():
    name = _cli("{id}_BC1:{1}_BC2:{2}_umi:rc({3})")
    assert name == "@x/1_BC1:AAAAA_BC2:GGGGG_umi:ATCACG"


def test_cli_empty_rename_rejected():
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            f1.write(f"@x/1\n{_R1}\n+\n{'I' * len(_R1)}\n")
            f2.write(f"@x/2\n{_R2}\n+\n{'I' * len(_R2)}\n")
        p = subprocess.run(
            [CUTSEQ, "-A", "INLINE", "-O", f"{td}/o", "--rename", "",
             f"{td}/t1.fq", f"{td}/t2.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode != 0
        assert "cannot be empty" in p.stderr
        assert "Traceback" not in p.stderr


def test_cli_bad_rename_clean_error():
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1:
            f1.write(f"@x/1\n{_R1}\n+\n{'I' * len(_R1)}\n")
        p = subprocess.run(
            [CUTSEQ, "-A", "INLINE", "-O", f"{td}/o", "--rename", "{id}_{9}",
             f"{td}/t1.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode != 0
        assert "cutseq:" in p.stderr        # our prefix, not a bare traceback
        assert "Traceback" not in p.stderr


def test_threaded_rename_deterministic():
    """-t N must not change output names (capture registry is safe across
    worker threads)."""
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            for k in range(300):
                f1.write(f"@x{k}/1\n{_R1}\n+\n{'I' * len(_R1)}\n")
                f2.write(f"@x{k}/2\n{_R2}\n+\n{'I' * len(_R2)}\n")
        tpl = "{id}_BC1:{1}_BC2:{2}_umi:rc({3})"

        def names(nthreads):
            p = subprocess.run(
                [CUTSEQ, "-A", "INLINE", "-O", f"{td}/o{nthreads}", "--rename", tpl,
                 "-t", str(nthreads), "-m", "10", f"{td}/t1.fq", f"{td}/t2.fq"],
                capture_output=True, text=True,
            )
            assert p.returncode == 0, p.stderr
            with gzip.open(f"{td}/o{nthreads}_trimmed_R1.fastq.gz", "rt") as f:
                out = []
                for i, line in enumerate(f):
                    if i % 4 == 0:
                        out.append(line.strip())
            return out

        a = names(1)
        b = names(4)
        assert a and a == b


# --- spatial (DBiT) interleaved 3' structure ---------------------------------
#
# m6A-ARTR-DBiT / Glori style, Read 1 (single-end):
#   handle : BCB linker2 BCA linker1 UMI polyT cDNA Tn5 R2H i5 P5
# (DBiT-seq spatial library; the m6A-ARTR-DBiT paper builds on this arm)
# where BCB/BCA are 8-nt spatial barcodes, UMI is 10 nt. P7/i7 sit 5' of the
# R1 primer and are NOT read. The 3' side must be processed in REVERSED
# written (physical) order, not grouped by token kind -- this was a real bug.


def _dbit_pair():
    import random

    random.seed(11)
    # Real DBiT-seq (Glori-style) oligo (Read 1 arm, from the protocol's
    # full double-stranded sequence). P7/i7 lie 5' of the Read-1 primer, so
    # they are NOT read; R1 starts at the PCR handle.
    handle = "CAAGCGTTGGCTTCTCGCATCT"          # PCR handle, 22
    linker2 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG" # 30
    linker1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT" # 30
    t5 = "CTGTCTCTTATACACATCT"                 # Tn5 mosaic end (read 2)
    r2h = "GACGCTGCCGACGA"                     # read-2 handle
    i5 = "TAGATCGC"
    p5 = "GTGTAGATCTCGGTGGTCGCCGTATCATT"
    bcB = "".join(random.choice("ACGT") for _ in range(8))
    bcA = "".join(random.choice("ACGT") for _ in range(8))
    umi = "".join(random.choice("ACGT") for _ in range(10))
    if umi[-1] == "T":
        umi = umi[:-1] + "A"
    cdna = "".join(random.choice("ACGT") for _ in range(30))
    # R1 5'->3': handle BCB linker2 BCA linker1 UMI polyT cDNA Tn5 R2H i5 P5
    r1 = handle + bcB + linker2 + bcA + linker1 + umi + "T" * 15 \
        + cdna + t5 + r2h + i5 + p5
    return r1, bcB, bcA, umi


def test_spatial_scheme_single_end_physical_order():
    from cutseq.grammar import _capture_registry

    r1, bcB, bcA, umi = _dbit_pair()
    # space-free scheme for Read 1 (single-end). P7/i7 are 5'-upstream of the
    # R1 primer (not read); BCB..polyT are the 5' side, the cDNA insert is
    # implicit between polyT and the Nextera Tn5/R2/i5/P5 3' arm (so ':' sits
    # right before the Tn5 arm).
    scheme = ("CAAGCGTTGGCTTCTCGCATCT"
              "N8ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
              "N8GTGGCCGATGTTTCGCATCGGCGTACGACT"
              "N10TTT...TTT"
              ":CTGTCTCTTATACACATCTGACGCTGCCGACGATAGATCGC"
              "GTGTAGATCTCGGTGGTCGCCGTATCATT")
    s = CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    mods[-1] = cs.renamer(paired=False, name_format="{id}_BCB:{1}_BCA:{2}_UMI:{3}")
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        read = SequenceRecord("x", r1, "I" * len(r1))
        info = ModificationInfo(read)
        for m in mods:
            out = m(read, info)
            if out is not None:
                read = out
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        _capture_registry.clear()
    assert read.name == f"x_BCB:{bcB}_BCA:{bcA}_UMI:{umi}"


# --- m6A-ARTR / DBiT spatial library (R2 = reverse complement of R1) ---------
#
# Barcode arm (handle + BCB + linker2 + BCA + linker1 + UMI) sits 5' of the
# insert and is the R1 read; R2 is the reverse complement (standard Illumina
# paired-end). The paired compile must walk the interleaved arm on R1's 5' in
# written order (adapters and captures interleaved), not kind-phase grouped.

_M6A_HANDLE = "CAAGCGTTGGCTTCTCGCATCT"
_M6A_L2 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
_M6A_L1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
_M6A_TSO = "AAGCAGTGGTATCAACGCAGAGT"


def _m6a_pair():
    import random

    random.seed(41)
    rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
    bcb, bca, umi = rnd(8), rnd(8), rnd(10)
    if umi[-1] == "T":
        umi = umi[:-1] + "A"
    cdna = rnd(30)
    r1 = _M6A_TSO + cdna
    r2 = _M6A_HANDLE + bcb + _M6A_L2 + bca + _M6A_L1 + umi + cdna
    return r1, r2, bcb, bca, umi


def test_paired_standard_colon_interleaved_barcode_arm():
    """Full ':' scheme must walk the interleaved handle->BCB->linker2->BCA->
    linker1->UMI arm on R1's 5' in written order (R2 = its reverse complement).
    """
    from cutseq.grammar import _capture_registry

    cdna_read, bc_read, bcb, bca, umi = _m6a_pair()
    # full scheme: barcode arm (R1, written order) : R2 primer (34 nt, reverse
    # complement -> cutseq trims the Tn5ME read-through CTGTCTCTTATACACATCT..)
    _R2_PRIMER = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    scheme = (_M6A_HANDLE + "N8" + _M6A_L2 + "N8" + _M6A_L1 + "N10"
              + "TTT...TTT" + ":" + _R2_PRIMER)
    s = CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=True, settings=s)
    mods[-1] = cs.renamer(paired=True, name_format="{id}_BCB:{1}_BCA:{2}_UMI:{3}")
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        read1 = SequenceRecord("b/1", bc_read, "I" * len(bc_read))   # R1: barcode arm
        read2 = SequenceRecord("c/2", cdna_read, "I" * len(cdna_read))  # R2: cDNA
        i1, i2 = ModificationInfo(read1), ModificationInfo(read2)
        for step in mods:
            if isinstance(step, tuple):
                m1, m2 = step
                n1 = m1(read1, i1) if m1 else read1
                n2 = m2(read2, i2) if m2 else read2
                read1, read2 = n1 or read1, n2 or read2
            else:
                step(read1, read2, i1, i2)
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        _capture_registry.clear()
    assert read1.name == f"b/1_BCB:{bcb}_BCA:{bca}_UMI:{umi}"
    # R1 barcode arm fully trimmed: only the 30 nt cDNA tail remains
    assert read1.sequence == bc_read[len(_M6A_HANDLE) + 8 + len(_M6A_L2)
                                + 8 + len(_M6A_L1) + 10:]


def test_dbitseq_builtin_renamed_from_m6aartr():
    """The spatial 50x50 barcode-arm scheme was misnamed 'M6AARTR'; it is
    DBiT-seq (deterministic barcoding in tissue, Cell 2020). The old name
    must remain a working alias pointing at the same scheme."""
    assert "DBITSEQ" in BUILDIN_ADAPTERS
    assert "M6AARTR" in BUILDIN_ADAPTERS          # kept as a legacy alias
    assert BUILDIN_ADAPTERS["M6AARTR"] == BUILDIN_ADAPTERS["DBITSEQ"]
    assert "N12AGTCGTACGCCGATGCGAAAC" in BUILDIN_ADAPTERS["DBITSEQ"]


def test_polytail_direction_auto_detected():
    """`B...B` needs no '^': its 5'/3' direction is inferred from the scheme
    layout. A run at the arm's read-start (first token, or right after the
    outer adapter) trims the read 5'; elsewhere (after captures/inner
    linkers) it trims the 3'. Leftmost-anchored on the 5' side."""
    from cutseq.grammar import (Poly5TailModifier, PolyTailModifier,
                                _Token, _mark_polytail_direction, tokenize)

    # dot-form parses as a plain polytail token, no marker needed
    assert [t.kind for t in tokenize("GGG...GGG")] == ["polytail"]

    # auto-detection of the trimmed end
    leading = [_Token("adp", "ACGTGTACA"), _Token("polytail", "G")]
    trailing = [_Token("adp", "ACGTGTACA"), _Token("capture", 4),
                _Token("polytail", "G")]
    _mark_polytail_direction(leading)
    _mark_polytail_direction(trailing)
    assert leading[1].options.get("five") is True      # leading -> 5'
    assert trailing[2].options.get("five") is None     # inner       -> 3'

    # the 5' modifier is leftmost-anchored (position 0)
    t5 = Poly5TailModifier("G", min_len=3)
    cases = [("GGGGGGACGT", "ACGT"),          # leading run trimmed
             ("ACGTACGGGGCGTT", "ACGTACGGGGCGTT"),  # internal: kept
             ("GGGGGGGGG", ""),               # all-G read
             ("GACGTACGT", "GACGTACGT")]      # single G (<min_len): kept
    for seq, want in cases:
        read = SequenceRecord("x", seq, "I" * len(seq))
        out = t5(read, ModificationInfo(read))
        assert out.sequence == want, (seq, out.sequence)

    # a trailing run still trims the 3' end
    t3 = PolyTailModifier("T", min_len=3)
    read = SequenceRecord("y", "ACGTACGTTTTTT", "I" * 13)
    assert t3(read, ModificationInfo(read)).sequence == "ACGTACG"


# --- 5' sequencing-primer args & auto-detection ------------------------------
def test_seq_primer_detection_db():
    from cutseq.primers import detect_5prime_primer

    hit = detect_5prime_primer("ACACTCTTTCCCTACACGACGCTCTTCCGATCTACGT")
    assert hit is not None and "TruSeq" in hit[1]
    # custom TSO / handle are not in the DB -> None (falls back to args)
    assert detect_5prime_primer("AAGCAGTGGTATCAACGCAGAGTAG") is None
    assert detect_5prime_primer("CAAGCGTTGGCTTCTCGCATCTACGT") is None


def test_seq_primer1_single_end_trims_tso():
    """--seq-primer1 trims R1's 5' TSO scaffold."""
    from cutseq.grammar import _capture_registry

    r1, r2, bcb, bca, umi = _m6a_pair()
    s = CutadaptConfig()
    s.seq_primer1 = _M6A_TSO
    cs = _build_scheme(":", s)   # `:` = keep the cDNA insert, no left tokens
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    mods[-1] = cs.renamer(paired=False, name_format="{id}")
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        read = SequenceRecord("x", r1, "I" * len(r1))
        info = ModificationInfo(read)
        for m in mods:
            out = m(read, info)
            if out is not None:
                read = out
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        _capture_registry.clear()
    assert read.sequence == r1[len(_M6A_TSO):]


def test_seq_primer2_single_end_trims_handle():
    """--seq-primer2 trims the 5' handle so the scheme can start at BCB."""
    from cutseq.grammar import _capture_registry

    r1, r2, bcb, bca, umi = _m6a_pair()
    s = CutadaptConfig()
    s.seq_primer2 = _M6A_HANDLE
    scheme = "N8" + _M6A_L2 + "N8" + _M6A_L1 + "N10"
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    mods[-1] = cs.renamer(paired=False, name_format="{id}_BCB:{1}_BCA:{2}_UMI:{3}")
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        read = SequenceRecord("x", r2, "I" * len(r2))
        info = ModificationInfo(read)
        for m in mods:
            out = m(read, info)
            if out is not None:
                read = out
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        _capture_registry.clear()
    assert read.name == f"x_BCB:{bcb}_BCA:{bca}_UMI:{umi}"
    assert read.sequence == r2[len(_M6A_HANDLE) + 8 + len(_M6A_L2)
                              + 8 + len(_M6A_L1) + 10:]


def test_primers_are_not_trimmed_from_reads():
    """Sequencing primers anneal upstream of each read and are NOT part of
    the read, so --r1-primer / --r2-primer must not inject any 5' trim."""
    from cutseq.run import _describe

    s = CutadaptConfig()
    s.r1_primer = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    s.r2_primer = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    cs = _build_scheme("ACGT:N8", s)
    assert all(t.kind != "adp" or t.value not in (s.r1_primer, s.r2_primer)
               for t in cs.left + cs.right)
    mods = _scheme_modifiers(cs, paired=True, settings=s)
    for m in mods:
        if isinstance(m, tuple):
            for mm, sub in ((m[0], "TCGTCGGCAGCGTC"),
                            (m[1], "GTCTCGTGGGCTCGG")):
                if mm:
                    assert sub not in _describe(mm)


def test_cli_rename_single_end():
    # single-end on the same scheme: capture2 & the inline are at R1's 3' end
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1:
            f1.write(f"@x\n{_R1}\n+\n{'I' * len(_R1)}\n")
        p = subprocess.run(
            [CUTSEQ, "-A", "INLINE", "-O", f"{td}/s", "--rename",
             "{id}_c1={1}_c2={2}_bc={3}", "-m", "10", f"{td}/t1.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        with gzip.open(f"{td}/s_trimmed_R1.fastq.gz", "rt") as f:
            name = next(f).strip()
    assert name.startswith("@x_c1=")


if __name__ == "__main__":
    import traceback

    failures = 0
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {nm}")
            except Exception:
                failures += 1
                print(f"FAIL {nm}")
                traceback.print_exc()
    sys.exit(1 if failures else 0)
