"""Tests for ``-A <file.yaml>`` scheme support (per-part settings).

Covers: YAML -> (orientation, left, right) parsing, the ``strand`` part acting
as the R1/R2 split, and per-part ``max_errors``/``min_overlap`` being applied to
the built modifiers independently on R1 vs R2.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")
R1 = str(ROOT / "test" / "input_R1.fq.gz")
R2 = str(ROOT / "test" / "input_R2.fq.gz")

sys.path.insert(0, str(ROOT))
from cutseq import grammar  # noqa: E402

PAIRED = """\
parts:
  - type: adapter
    seq: AGTTCTACAGTCCGACGATC
    max_errors: 0.2
    min_overlap: 10
  - type: insert
    value: "+"
  - type: adapter
    seq: AGATCGGAAGAGCACACGTC
    max_errors: 0.1
    min_overlap: 3
"""

SINGLE = """\
parts:
  - type: adapter
    seq: ACACGACGCTCTTCCGATCT
    max_errors: 0.2
  - type: mask
    length: 3
  - type: capture
    length: 8
    label: umi
  - type: inline
    seq: atcacg
    label: bc1
    max_errors: 1
"""


def _write(tmp, name, text):
    p = Path(tmp) / name
    p.write_text(text)
    return str(p)


def test_parse_paired_insert_part():
    import yaml

    data = yaml.safe_load(PAIRED)
    orientation, left, right = grammar.parse_scheme_parts(data)
    assert orientation == "+"
    assert [(t.kind, t.value) for t in left] == [("adp", "AGTTCTACAGTCCGACGATC")]
    assert [(t.kind, t.value) for t in right] == [("adp", "AGATCGGAAGAGCACACGTC")]


def test_parse_single_no_strand():
    import yaml

    data = yaml.safe_load(SINGLE)
    orientation, left, right = grammar.parse_scheme_parts(data)
    assert orientation is None
    assert right == []
    assert [(t.kind, t.value) for t in left] == [
        ("adp", "ACACGACGCTCTTCCGATCT"),
        ("mask", 3),
        ("capture", 8),
        ("inline", "atcacg"),
    ]
    assert left[2].label == "umi"
    assert left[3].label == "bc1"


def test_per_part_errors_apply_to_each_side():
    import yaml

    data = yaml.safe_load(PAIRED)
    orientation, left, right = grammar.parse_scheme_parts(data)
    m1, m2, _ = grammar.build_modifiers_from_parts(
        orientation, left, right, paired=True
    )
    # read-plan routing: R1's own-arm adapter is the left adapter, and its
    # read-through gate is the right (p7) adapter both honour their per-part
    # options; the p7 site is upstream/absent from R2 (it is R2's boundary).
    steps = m1[0].step_adapters()
    agtt = [a for (s, a) in steps
            if getattr(a, "sequence", None) == "AGTTCTACAGTCCGACGATC"][0]
    assert agtt.max_error_rate == 0.2
    assert agtt.min_overlap >= 10  # requested floor (cutadapt clamps upward)
    p7 = [a for (s, a) in steps
          if getattr(a, "sequence", None) == "AGATCGGAAGAGCACACGTC"][0]
    assert p7.max_error_rate == 0.1
    assert p7.min_overlap >= 3


def test_strand_marker_rcs_r2():
    # The read-2 mechanism reverse-complements the written top-strand adapter
    # before matching it, so R2's read-through of the LEFT arm is the rc of the
    # left adapter (and the right outer site is R2's upstream boundary, not
    # matched on R2 itself).
    import yaml

    data = yaml.safe_load(PAIRED)
    orientation, left, right = grammar.parse_scheme_parts(data)
    m1, m2, _ = grammar.build_modifiers_from_parts(
        orientation, left, right, paired=True
    )
    r2_seqs = [getattr(a, "sequence", None) for (s, a) in m2[0].step_adapters()]
    # The read-2 mechanism rc's the written top-strand left adapter before it is
    # matched as R2's read-through gate.
    assert grammar._rc("AGTTCTACAGTCCGACGATC") in r2_seqs


def test_cli_accepts_yaml_file_paired():
    with tempfile.TemporaryDirectory() as td:
        y = _write(td, "scheme.yaml", PAIRED)
        p = subprocess.run(
            [CUTSEQ, "-A", y, "-n", R1, R2], capture_output=True, text=True
        )
        assert p.returncode == 0, p.stderr
        # read-plan routing: each read is an arm-local executor, visible in -n.
        assert "PlanExecutor" in p.stderr + p.stdout


def test_cli_accepts_yaml_file_single():
    with tempfile.TemporaryDirectory() as td:
        y = _write(td, "scheme.yaml", SINGLE)
        p = subprocess.run(
            [CUTSEQ, "-A", y, "-n", R1], capture_output=True, text=True
        )
        assert p.returncode == 0, p.stderr
        assert "ConditionalCutter(length=8" in p.stderr + p.stdout
        # The YAML 'label: umi' capture should be visible in dry-run output.
        assert "name=umi" in p.stderr + p.stdout


def test_cli_runs_paired_yaml():
    with tempfile.TemporaryDirectory() as td:
        y = _write(td, "scheme.yaml", PAIRED)
        prefix = str(Path(td) / "out")
        p = subprocess.run(
            [CUTSEQ, "-A", y, "-O", prefix, R1, R2],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        assert (Path(td) / "out_trimmed_R1.fastq.gz").exists()


def test_full_adapter_settings_flow():
    """All cutadapt adapter kwargs (indels, wildcards, force_anywhere) pass through."""
    import yaml

    data = yaml.safe_load("""\
parts:
  - type: adapter
    seq: AGTTCTACAGTCCGACGATC
    max_errors: 0.2
    min_overlap: 10
    indels: false
    read_wildcards: true
  - type: insert
    value: "+"
  - type: adapter
    seq: AGATCGGAAGAGCACACGTC
    force_anywhere: true
    adapter_wildcards: false
""")
    orientation, left, right = grammar.parse_scheme_parts(data)
    m1, m2, _ = grammar.build_modifiers_from_parts(
        orientation, left, right, paired=True
    )
    # read-plan: R1's own-arm left adapter (indels/wildcards/max_errors) and its
    # read-through right adapter (force_anywhere / adapter_wildcards) honour the
    # per-part options on the exact step adapters.
    steps = m1[0].step_adapters()
    a1 = [a for (s, a) in steps
          if getattr(a, "sequence", None) == "AGTTCTACAGTCCGACGATC"][0]
    assert a1.max_error_rate == 0.2
    assert a1.min_overlap >= 10  # requested floor (cutadapt clamps upward)
    assert a1.indels is False
    assert a1.read_wildcards is True
    a2 = [a for (s, a) in steps
          if getattr(a, "sequence", None) == "AGATCGGAAGAGCACACGTC"][0]
    assert a2.adapter_wildcards is False
    assert a2._force_anywhere is True


def test_bare_list_yaml_form():
    """A scheme file may be a bare list (no top-level 'parts:' wrapper)."""
    import yaml

    data = yaml.safe_load("""\
- type: adapter
  seq: AGTTCTACAGTCCGACGATC
- type: insert
  value: "+"
- type: adapter
  seq: AGATCGGAAGAGCACACGTC
""")
    orientation, left, right = grammar.parse_scheme_parts(data)
    assert orientation == "+"
    assert [(t.kind, t.value) for t in left] == [("adp", "AGTTCTACAGTCCGACGATC")]
    assert [(t.kind, t.value) for t in right] == [("adp", "AGATCGGAAGAGCACACGTC")]


def test_unknown_part_type_errors():
    import yaml

    data = yaml.safe_load("parts:\n  - type: bogus\n    seq: ACGT\n")
    try:
        grammar.parse_scheme_parts(data)
    except ValueError as e:
        assert "unknown part type" in str(e)
    else:
        raise AssertionError("expected ValueError")


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
