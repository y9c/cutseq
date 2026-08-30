"""Tests for CLI flags reconnected to the grammar pipeline.

These cover the behaviour that the grammar refactor originally dropped:
  - --min-quality  -> a QualityTrimmer is present in the modifier chain
  - --with-rname-suffix -> SuffixRemover steps are prepended
  - --auto-rc      -> ReverseComplementModifier added for '-' schemes (single-end)
  - --conditional-cutter / --force-trim-min-length -> capture/mask only cut
    when an adapter matched or the read is long enough
  - metrics: quality trimming should route short reads to the short file
"""

import gzip
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")
R1 = str(ROOT / "test" / "input_R1.fq.gz")
R2 = str(ROOT / "test" / "input_R2.fq.gz")

sys.path.insert(0, str(ROOT))
from cutseq import grammar  # noqa: E402


def _dry(scheme, *extra, single=True, use_r2=False):
    args = [CUTSEQ, "-A", scheme, R1]
    if use_r2:
        args.append(R2)
    args += ["-n", *extra]
    out = subprocess.run(args, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stderr + "\n" + out.stdout


def _modsteps(output):
    """Return the lines of the dry-run step list."""
    return [line for line in output.splitlines() if line.startswith("Step")]


def test_min_quality_adds_quality_trimmer():
    out = _dry("SMALLRNA", "-q", "30", use_r2=True)
    steps = "\n".join(_modsteps(out))
    assert "QualityTrimmer" in steps
    assert "cutoff_back=30" in steps


def test_quality_trimming_sends_short_reads_to_short_file():
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        prefix = str(Path(td) / "out")
        p = subprocess.run(
            [CUTSEQ, "-A", "SMALLRNA", "-q", "40", "-O", prefix, R1, R2],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        discard = Path(td) / "out_discard_R1.fastq.gz"
        trimmed = Path(td) / "out_trimmed_R1.fastq.gz"
        assert discard.exists() and trimmed.exists()
        with gzip.open(discard, "rt") as f:
            names = [next(f) for _ in range(4)]
        with gzip.open(discard, "rt") as f:
            discard_reads = sum(1 for _ in f) // 4
        with gzip.open(trimmed, "rt") as f:
            trimmed_reads = sum(1 for _ in f) // 4
        assert discard_reads > 0, "quality trimming should have discarded some reads"
        assert discard_reads + trimmed_reads == 10000
        assert any("reason=too_short" in n for n in names)


def test_rname_suffix_prepends_suffix_remover():
    out = _dry("SMALLRNA", "--with-rname-suffix", use_r2=True)
    steps = _modsteps(out)
    assert steps and "SuffixRemover" in steps[0], steps


def test_auto_rc_single_antisense_adds_rc_modifier():
    # '-' is the antisense strand marker: for single-end it triggers RC.
    out = _dry("TAKARAV2", "--auto-rc")
    steps = "\n".join(_modsteps(out))
    assert "ReverseComplement" in steps


def test_auto_rc_single_nonminus_warns():
    out = _dry("UNSTRANDED", "--auto-rc")
    assert "ignoring --auto-rc" in out


def test_auto_rc_paired_warns():
    out = _dry("SMALLRNA", "--auto-rc", use_r2=True)
    assert "--auto-rc is ignored for paired-end" in out


def test_conditional_cutter_keeps_short_unmatched_capture():
    # With --conditional-cutter and a small force-trim-min-length, the
    # 3'-side capture cutter carries conditional=True / force_trim_min_length
    # so a capture on a short, adapter-free read is NOT cut.
    mods, _, orient = grammar.build_modifiers(
        "ACGT+N8", conditional_cutter=True, force_trim_min_length=50
    )
    assert orient == "+"
    cap = mods[1]
    assert cap.conditional is True
    assert cap.force_trim_min_length == 50


def test_conditional_cutter_disabled_cuts_always():
    mods, _, _ = grammar.build_modifiers(
        "ACGT+N8", conditional_cutter=False
    )
    assert mods[1].conditional is False


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
