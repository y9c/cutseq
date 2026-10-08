"""Tests for the extensible engine API: CompiledScheme, the token-kind
registry, and read-name customization (--name-format / --rename)."""

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cutseq"))

from cutseq import grammar  # noqa: E402


def _rc(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def test_compiled_scheme_has_captures():
    cs = grammar.CompiledScheme("+", [], [grammar._Token("capture", 8)])
    assert cs.has_captures is True
    cs2 = grammar.CompiledScheme("+", [], [grammar._Token("mask", 3)])
    assert cs2.has_captures is False


def test_compiled_scheme_modifiers_single():
    """Structural/API contract of ``CompiledScheme.modifiers(paired=False)``: a
    single-end scheme yields a 5' front-adapter modifier plus a 3' UMI-capture
    modifier (and no R2 column).  This is introspection of the compiled chain,
    not a run of the pipeline."""
    cs = grammar.CompiledScheme(None, grammar.parse_scheme("ACGT+N8")[1],
                                grammar.parse_scheme("ACGT+N8")[2])
    mods, r2 = cs.modifiers(paired=False)
    assert r2 == []
    # adapter (front) + 3' UMI capture
    assert len(mods) == 2
    assert mods[1].capture is True


def test_compiled_scheme_modifiers_paired():
    """A primer-boundary scheme (ECLIP10-style) compiles to per-read plans that
    actually trim the UMI/masks and preserve the insert.

    The outer p5/p7 adapters are genuine sequencing-site binders (upstream,
    absent from the reads), so each read starts inside them.  Verified by
    running the compiled R1/R2 modifiers on a hand-built read pair and checking
    the trimmed sequences and the captured UMI literal — NOT by counting
    internal executor steps (which is an implementation detail, not behavior).
    """
    import random

    from dnaio import SequenceRecord
    from cutadapt.info import ModificationInfo

    _, left, right = grammar.parse_scheme("ACACGACGCTCTTCCGATCTXX-XNNNNNNNNNNAGATCGGAAGAGCACACGTC")
    cs = grammar.CompiledScheme("-", left, right)
    m1, m2 = cs.modifiers(paired=True)
    ren = cs.renamer(paired=True, name_format="{id}_r1:{r1.1}_r2:{r2.1}")
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        rnd = random.Random(4)
        rnd_nt = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
        x2, ins, x1, n10 = rnd_nt(2), rnd_nt(30), rnd_nt(1), rnd_nt(10)
        r1p = "ACACGACGCTCTTCCGATCT"
        p7 = "AGATCGGAAGAGCACACGTC"
        # R1 starts after the p5 (r1p) site: mask + insert + mask + UMI + p7
        # read-through; R2 starts after the p7 site on the bottom strand.
        r1_seq = x2 + ins + x1 + n10 + p7
        r2_seq = _rc(n10) + _rc(x1) + _rc(ins) + _rc(x2) + _rc(r1p)
        r1 = SequenceRecord("x/1", r1_seq, "I" * len(r1_seq))
        r2 = SequenceRecord("x/2", r2_seq, "I" * len(r2_seq))
        i1, i2 = ModificationInfo(r1), ModificationInfo(r2)
        for mm in m1:
            out = mm(r1, i1)
            r1 = out or r1
        for mm in m2:
            out = mm(r2, i2)
            r2 = out or r2
        ren(r1, r2, i1, i2)
        assert r1.sequence == ins
        assert r2.sequence == _rc(ins)
        assert r1.name == f"x/1_r1:{n10}_r2:{_rc(n10)}"
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def test_make_renamer_custom_format():
    r = grammar.make_renamer(False, has_captures=True, name_format="{id}|{cut_suffix}")
    assert r._template == "{id}|{cut_suffix}"


def test_make_renamer_name_format_template():
    # the --rename template is stored verbatim for introspection
    orientation, left, right = grammar.parse_scheme("ACGT+N8")
    cs = grammar.CompiledScheme(orientation, left, right)
    r = cs.renamer(paired=True, name_format="{id}:{1}")
    assert r._template == "{id}:{1}"


def test_make_renamer_defaults_match_legacy():
    assert grammar.make_renamer(True, has_captures=True)._template == \
        "{id}_{r1.cut_prefix}{r2.cut_prefix}"
    assert grammar.make_renamer(True, has_captures=False)._template == "{id}"
    assert grammar.make_renamer(False, has_captures=True)._template == \
        "{id}_{cut_prefix}{cut_suffix}"


def test_emitter_registry_covers_all_kinds():
    for kind in ("adp", "back", "inline", "capture", "mask", "polytail"):
        assert kind in grammar._EMITTERS
        five, three, se_three = grammar._EMITTERS[kind]
        assert five is not None
    # every registry kind appears in the single-end + paired phase tables
    assert set(grammar._EMITTERS) == set(grammar._SINGLE_KINDS)
    assert {k for k, _ in grammar._PAIRED_PHASES} == set(grammar._EMITTERS)


def test_cli_name_format_flag():
    import subprocess
    cutseq = str(ROOT / ".venv" / "bin" / "cutseq")
    r = subprocess.run([cutseq, "--help"], capture_output=True, text=True)
    assert "--name-format" in r.stdout
    assert "--rename" in r.stdout
    assert "-R" in r.stdout


def test_cli_short_R_rename_flag():
    import subprocess
    import gzip
    cutseq = str(ROOT / ".venv" / "bin" / "cutseq")
    r1 = str(ROOT / "test" / "input_R1.fq.gz")
    with tempfile.TemporaryDirectory() as td:
        p = subprocess.run(
            [cutseq, "-A", "ECLIP10", "-R", "{id}|{1}",
             "-O", f"{td}/cam", r1],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        with gzip.open(f"{td}/cam_trimmed_R1.fastq.gz", "rt") as fh:
            first = next(fh).strip()
        assert "|" in first


def test_cli_name_format_end_to_end():
    import subprocess
    import gzip
    cutseq = str(ROOT / ".venv" / "bin" / "cutseq")
    r1 = str(ROOT / "test" / "input_R1.fq.gz")
    with tempfile.TemporaryDirectory() as td:
        p = subprocess.run(
            [cutseq, "-A", "ECLIP10", "--name-format", "{id}|{cut_suffix}",
             "-O", f"{td}/cam", r1],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        out = f"{td}/cam_trimmed_R1.fastq.gz"
        with gzip.open(out, "rt") as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    assert line.startswith("@"), line
                    assert "|" in line
                    break


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
