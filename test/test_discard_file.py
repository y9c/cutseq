"""Tests for the unified --discard-file output with reason-tagged read names.

Covers:
  - reason=too_short  (shorter than --min-length after trimming)
  - reason=too_many_n (exceeds --max-n)
  - reason=low_quality (mean Phred quality below --min-avg-quality)
  - reason=no_barcode (missing an expected inline barcode, --ensure-inline-barcode)
  - reads with a valid inline barcode pass through to the trimmed output
  - single-end and paired-end pipelines both write the discard file
  - json_report emits discard1/discard2 instead of short/untrimmed
"""

import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")
R1 = str(ROOT / "test" / "input_R1.fq.gz")
R2 = str(ROOT / "test" / "input_R2.fq.gz")

sys.path.insert(0, str(ROOT))


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_pair(td, reads):
    """Write a 4-record FASTQ pair. Each entry: (name, r1_seq, r2_seq)."""
    r1p = Path(td) / "t1.fq"
    r2p = Path(td) / "t2.fq"
    with open(r1p, "w") as f1, open(r2p, "w") as f2:
        for name, s1, s2 in reads:
            f1.write(f"@{name}/1\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@{name}/2\n{s2}\n+\n{'I' * len(s2)}\n")
    return str(r1p), str(r2p)


def _names(fastq_gz):
    with gzip.open(fastq_gz, "rt") as f:
        return [line.strip() for line in f if line.startswith("@")]


def _build_reads(prefix, insert, bc_phys="ATCACG", p5=None, p7=None):
    """Build correctly-oriented R1/R2 for scheme P5 + insert + bc + P7."""
    p5 = p5 or "AGTTCTACAGTCCGACGATC"
    p7 = p7 or "AGATCGGAAGAGCACACGTC"
    bcrc = _rc(bc_phys)
    r1 = p5 + insert + bc_phys + p7
    r2 = _rc(p7) + bcrc + _rc(insert) + _rc(p5)
    return (prefix, r1, r2)


def _inline_scheme(bc_phys="ATCACG"):
    """Custom grammar scheme with an inline barcode between the two primers."""
    return ("AGTTCTACAGTCCGACGATC+" + bc_phys.lower()
            + "AGATCGGAAGAGCACACGTC")


def _run_paired(td, scheme, reads, extra=()):
    r1, r2 = _write_pair(td, reads)
    prefix = str(Path(td) / "out")
    p = subprocess.run(
        [CUTSEQ, "-A", scheme, "-O", prefix, "-m", "50", *extra, r1, r2],
        capture_output=True, text=True,
    )
    assert p.returncode == 0, p.stderr
    return p


def test_discard_file_tags_all_reasons_paired():
    with tempfile.TemporaryDirectory() as td:
        ins = "ACGT" * 25
        reads = [
            _build_reads("ok", ins),
            _build_reads("no_bc", ins, bc_phys="TTTTTT"),
            _build_reads("short", "ACG"),
            _build_reads("manyn", "NNNNN" * 5),
            _build_reads("lowq", ins),
        ]
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            for name, s1, s2 in reads:
                q = "!" if name == "lowq" else "I"
                f1.write(f"@{name}/1\n{s1}\n+\n{q * len(s1)}\n")
                f2.write(f"@{name}/2\n{s2}\n+\n{q * len(s2)}\n")
        prefix = str(Path(td) / "out")
        p = subprocess.run(
            [CUTSEQ, "-A", _inline_scheme(), "-O", prefix, "-m", "50",
             "--ensure-inline-barcode", "--max-n", "5",
             "--min-avg-quality", "30", "-q", "0",
             f"{td}/t1.fq", f"{td}/t2.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        discard = Path(td) / "out_discard_R1.fastq.gz"
        trimmed = Path(td) / "out_trimmed_R1.fastq.gz"
        assert discard.exists() and trimmed.exists()
        names = _names(str(discard))
        assert any("ok/1" in n for n in names) is False
        assert any("reason=too_short" in n and "short/1" in n for n in names)
        assert any("reason=too_many_n" in n and "manyn/1" in n for n in names)
        assert any("reason=low_quality" in n and "lowq/1" in n for n in names)
        assert any("reason=no_barcode" in n and "no_bc/1" in n for n in names)
        trimmed_names = _names(str(trimmed))
        assert any("ok/1" in n for n in trimmed_names)


def test_discard_requires_ensure_inline_barcode_for_untrimmed():
    with tempfile.TemporaryDirectory() as td:
        ins = "ACGT" * 25
        reads = [
            _build_reads("ok", ins),
            _build_reads("no_bc", ins, bc_phys="TTTTTT"),
        ]
        # Without --ensure-inline-barcode the no-barcode reads must NOT be
        # tagged; they go to the trimmed output (no discard for barcode).
        _run_paired(td, _inline_scheme(), reads)
        discard = Path(td) / "out_discard_R1.fastq.gz"
        names = _names(str(discard)) if discard.exists() else []
        assert not any("reason=no_barcode" in n for n in names)


def test_discard_file_single_end():
    with tempfile.TemporaryDirectory() as td:
        # single-end: just one file, min-length discards short reads
        ins = "ACGT" * 25
        p5 = "AGTTCTACAGTCCGACGATC"
        p7 = "AGATCGGAAGAGCACACGTC"
        reads = [
            ("ok", p5 + ins + p7, ""),
            ("short", p5 + "ACG" + p7, ""),
        ]
        r1p, _ = _write_pair(td, reads)
        prefix = str(Path(td) / "out")
        p = subprocess.run(
            [CUTSEQ, "-A", "SMALLRNA", "-O", prefix, "-m", "50", r1p],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        discard = Path(td) / "out_discard_R1.fastq.gz"
        assert discard.exists()
        names = _names(str(discard))
        assert any("reason=too_short" in n and "short" in n for n in names)


def test_json_report_has_discard_fields():
    with tempfile.TemporaryDirectory() as td:
        ins = "ACGT" * 25
        reads = [_build_reads("ok", ins), _build_reads("short", "ACG")]
        r1, r2 = _write_pair(td, reads)
        prefix = str(Path(td) / "out")
        jsonf = str(Path(td) / "report.json")
        p = subprocess.run(
            [CUTSEQ, "-A", "SMALLRNA", "-O", prefix, "-m", "50",
             "--json-file", jsonf, r1, r2],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr
        import json

        with open(jsonf) as f:
            data = json.load(f)
        out = data["output"]
        assert "discard1" in out and "discard2" in out
        assert "short1" not in out and "untrimmed1" not in out


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
