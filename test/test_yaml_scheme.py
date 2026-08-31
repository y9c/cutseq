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
    assert m1[0].adapters[0].max_error_rate == 0.2
    assert m1[0].adapters[0].min_overlap == 10
    assert m2[0].adapters[0].max_error_rate == 0.1
    assert m2[0].adapters[0].min_overlap == 3


def test_strand_marker_rcs_r2():
    # The R2 adapter must be reverse-complemented before matching.
    import yaml

    data = yaml.safe_load(PAIRED)
    orientation, left, right = grammar.parse_scheme_parts(data)
    m1, m2, _ = grammar.build_modifiers_from_parts(
        orientation, left, right, paired=True
    )
    assert m2[0].adapters[0].sequence == "GACGTGTGCTCTTCCGATCT"  # RC of AGATCGGAAGAGCACACGTC


def test_cli_accepts_yaml_file_paired():
    with tempfile.TemporaryDirectory() as td:
        y = _write(td, "scheme.yaml", PAIRED)
        p = subprocess.run(
            [CUTSEQ, "-A", y, "-n", R1, R2], capture_output=True, text=True
        )
        assert p.returncode == 0, p.stderr
        assert "AdapterCutter" in p.stderr + p.stdout


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
    a1, a2 = m1[0].adapters[0], m2[0].adapters[0]
    assert a1.max_error_rate == 0.2
    assert a1.min_overlap == 10
    assert a1.indels is False
    assert a1.read_wildcards is True
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
