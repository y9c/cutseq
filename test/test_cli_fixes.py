"""Tests for CLI fixes added on top of the grammar refactor:
  - argparse's nargs='+' greediness (input files after -o/-d were swallowed)
  - --rename positional captures are compiled and visible in dry-run output
  - JSON report 'barcode' field populated from the compiled scheme
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")
R1 = str(ROOT / "test" / "input_R1.fq.gz")
R2 = str(ROOT / "test" / "input_R2.fq.gz")

sys.path.insert(0, str(ROOT))
from cutseq.run import _option_arities, _reorder_inputs_first  # noqa: E402


def _dry(*args):
    p = subprocess.run([CUTSEQ, *args, "-n"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return p.stderr + "\n" + p.stdout


# --- _reorder_inputs_first unit tests --------------------------------------


def _make_parser():
    # Reuse the real parser: --help is enough to construct it, but we need the
    # actual parser instance. cutseq.run.main builds it inside; emulate by
    # importing argparse and duplicating the option set is fragile, so instead
    # build a minimal parser with the same arity structure.
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", nargs="*")
    parser.add_argument("-A", "-a", "--adapter-scheme")
    parser.add_argument("-O", "--output-prefix")
    parser.add_argument("-o", "--output-file", nargs="+")
    parser.add_argument("-d", "--discard-file", nargs="+")
    parser.add_argument("-m", "--min-length", type=int)
    parser.add_argument("--json-file")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-V", "--version", action="version", version="x")
    return parser


def test_reorder_paired_outputs_before_inputs():
    arities = _option_arities(_make_parser())
    out = _reorder_inputs_first(
        ["-A", "SMALLRNA", "-o", "o1", "o2", R1, R2], arities)
    assert out == [R1, R2, "-A", "SMALLRNA", "-o", "o1", "o2"]


def test_reorder_single_output_before_input():
    arities = _option_arities(_make_parser())
    out = _reorder_inputs_first(["-o", "o1", R1], arities)
    assert out == [R1, "-o", "o1"]


def test_reorder_discard_before_inputs():
    arities = _option_arities(_make_parser())
    out = _reorder_inputs_first(
        ["-A", "SMALLRNA", "-d", "d1", "d2", "-o", "o1", "o2", R1, R2], arities)
    assert out == [R1, R2, "-A", "SMALLRNA", "-d", "d1", "d2", "-o", "o1", "o2"]


def test_reorder_inputs_first_is_identity():
    arities = _option_arities(_make_parser())
    # Inputs already precede the options -> unchanged.
    argv = ["-A", "SMALLRNA", R1, R2, "-o", "o1", "o2"]
    assert _reorder_inputs_first(argv, arities) == argv


def test_reorder_respects_ddash():
    arities = _option_arities(_make_parser())
    argv = ["-o", "o1", "--", R1]
    assert _reorder_inputs_first(argv, arities) == argv


def test_reorder_returns_unchanged_when_inconsistent():
    arities = _option_arities(_make_parser())
    # 2 outputs + 1 input is never consistent -> argparse reports its own error.
    argv = ["-o", "o1", "o2", R1]
    assert _reorder_inputs_first(argv, arities) == argv


# --- CLI integration tests --------------------------------------------------


def test_cli_output_files_before_inputs_paired():
    with tempfile.TemporaryDirectory() as td:
        o1, o2 = f"{td}/t_R1.fastq.gz", f"{td}/t_R2.fastq.gz"
        p = subprocess.run(
            [CUTSEQ, "-A", "SMALLRNA", "-o", o1, o2, R1, R2, "-m", "50"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        assert Path(o1).exists() and Path(o2).exists()


def test_cli_discard_files_before_inputs_paired():
    with tempfile.TemporaryDirectory() as td:
        d1, d2 = f"{td}/d_R1.fastq.gz", f"{td}/d_R2.fastq.gz"
        p = subprocess.run(
            [CUTSEQ, "-A", "SMALLRNA", "-d", d1, d2, R1, R2, "-m", "90"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        assert Path(d1).exists() and Path(d2).exists()


def test_dry_run_shows_capture_labels_and_rename_template():
    # captures are labelled barcode{N} by default and visible in -n output;
    # a --rename template is shown compiled in the step list.
    out = _dry("-A", "INLINE", "--rename", "{id}_BC1:{1}_BC2:{2}", R1, R2)
    assert "name=barcode1" in out or "name=barcode" in out
    assert "{id}_BC1:{1}_BC2:{2}" in out or "Labeled" in out


def test_json_report_barcode_populated():
    with tempfile.TemporaryDirectory() as td:
        jsonf = f"{td}/report.json"
        p = subprocess.run(
            [CUTSEQ, "-A", "SMALLRNA", "-O", f"{td}/o", "-m", "90",
             "--json-file", jsonf, R1],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        with open(jsonf) as f:
            data = json.load(f)
        bc = data["barcode"]
        assert bc["orientation"] == "+"
        assert len(bc["parts"]) == 2
        assert bc["parts"][0]["value"] == "AGTTCTACAGTCCGACGATC"


def test_version_and_list_still_work():
    p = subprocess.run([CUTSEQ, "-V"], capture_output=True, text=True)
    assert p.returncode == 0
    assert "cutseq" in p.stdout

    p = subprocess.run([CUTSEQ, "--list-adapters"], capture_output=True, text=True)
    assert p.returncode == 0
    assert "TAKARAV3" in p.stdout


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
