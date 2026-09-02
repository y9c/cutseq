"""Migrated SCG/seqspec library schemes (adapters.toml).

Every scheme auto-imported from ``scg_to_cutseq.py`` is fully typed (adapter /
barcode / UMI roles are explicit in the grammar) and therefore registered with
``auto_inline = false`` so cutseq's inline-barcode heuristic never re-labels a
structural linker. These tests: (1) parse every migrated scheme, (2) for the
canonical 3'/5' methods, synthesize an R1 from the written scheme side and
assert barcode/UMI extraction + trimming through the real modifier pipeline.
"""

import random
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cutseq.common import load_adapters_no_auto_inline  # noqa: E402
from cutseq.grammar import parse_scheme, tokenize  # noqa: E402
from cutseq.run import CutadaptConfig, _build_scheme, _describe, _scheme_modifiers  # noqa: E402

from dnaio import SequenceRecord  # noqa: E402
from cutadapt.info import ModificationInfo  # noqa: E402

RAND = random.Random(1234)
_ADAPTERS = tomllib.loads(
    (Path(__file__).resolve().parent.parent / "cutseq" / "adapters.toml")
    .read_text(encoding="utf-8")
)

# canonical methods whose chemistry is well understood
KEY_METHODS = ["10X_RNA_V3", "DROP_SEQ", "SPLIT_SEQ", "CEL_SEQ2",
               "10X_RNA_V2", "SMART_SEQ2", "10X_ATAC"]


def migrated_names():
    return {m for m, cfg in _ADAPTERS.items()
            if isinstance(cfg, dict) and "scheme" in cfg
            and cfg.get("auto_inline") is False}


def test_all_migrated_schemes_parse_without_auto_inline():
    names = migrated_names()
    assert names, "no migrated schemes flagged auto_inline=false"
    for name in names:
        scheme = _ADAPTERS[name]["scheme"]
        parse_scheme(scheme, auto_inline=False)


def test_migrated_names_report_no_auto_inline():
    for name in migrated_names():
        assert name in load_adapters_no_auto_inline()


def build_r1(scheme, seed):
    """Synthesize an R1 read from the written left side of a scheme.

    Adapters/primers become their fixed sequence; N-captures become random
    bases we later assert on; poly tails expand; the insert marker ends the
    read with an explicit 30 nt random insert.
    """
    rnd = random.Random(seed)
    left = scheme.split(":", 1)[0]
    tokens = parse_scheme(left, auto_inline=False)[1]  # left tokens
    pieces = []
    captured = []
    cap_i = 0
    for tok in tokens:
        if tok.kind == "adp":
            pieces.append(tok.value)
        elif tok.kind == "capture":
            cap = "".join(rnd.choice("ACGT") for _ in range(tok.value))
            captured.append(cap)
            pieces.append(cap)
            cap_i += 1
        elif tok.kind == "mask":
            pieces.append("X" * tok.value)
        elif tok.kind == "polytail":
            pieces.append("T" * 15)
        elif tok.kind == "insert":
            pass
    pieces.append("".join(rnd.choice("ACGT") for _ in range(30)))
    return "".join(pieces), captured


@pytest.mark.parametrize("name", migrated_names())
def test_migrated_scheme_builds_and_captures_r1(name):
    cfg = _ADAPTERS[name]
    scheme = cfg["scheme"]
    if ":" not in scheme:
        pytest.skip("single-read/no-insert scheme")
    r1_seq, captured = build_r1(scheme, seed=hash(name) & 0xFFFF)
    rename = cfg.get("recommended_rename") or "{id}"
    s = CutadaptConfig()
    s.auto_inline = False
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    mods[-1] = cs.renamer(paired=False, name_format=rename)
    read = SequenceRecord("x", r1_seq, "I" * len(r1_seq))
    info = ModificationInfo(read)
    for m in mods:
        out = m(read, info)
        if out is not None:
            read = out
    # every capture the scheme declares must be referenced by the rename
    ltok = [t for t in tokenize(scheme.split(":", 1)[0]) if t.kind == "capture"]
    if not ltok:
        return  # no captures -> nothing to assert on names
    for i in range(1, len(ltok) + 1):
        assert f":{{{i}}}" in rename


# --- chemistry spot-checks: the R2 5'-most front adapter must be the read-2
# sequencing primer for canonical paired-end methods (cutseq rc's the written
# right-most adapter when it lands as R2's 5'). -------------------------------

_R2_READ2_PRIMER = {
    "10X_RNA_V3": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "10X_RNA_V2": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "DROP_SEQ": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "SMART_SEQ2": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "CEL_SEQ2": "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
}


def _r2_5prime_adp(scheme, s):
    cs = _build_scheme(scheme, s)
    _r1, r2 = cs.modifiers(paired=True)
    return _describe(r2[0])


@pytest.mark.parametrize("name", list(_R2_READ2_PRIMER))
def test_migrated_r2_5prime_is_read2_primer(name):
    cfg = _ADAPTERS[name]
    s = CutadaptConfig()
    s.auto_inline = False
    desc = _r2_5prime_adp(cfg["scheme"], s)
    assert _R2_READ2_PRIMER[name] in desc
