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
