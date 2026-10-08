#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Physically-annotated behavioral tests for the corrected built-in example
schemes in ``cutseq/adapters.toml``.

These are source-of-truth corrections for the previously-broken built-ins
(SMART_SEQ3, 10X_RNA_ATAC, 10X_RNA_ATAC_MULTI, SHARE_SEQ).  Each test builds a
hand-constructed read pair whose molecular geometry is known *independently* of
cutseq's resolver: the top-strand 5'->3' molecule is drawn from the corrected
scheme, R1 is read off the top strand and R2 is read off the bottom strand (so
R2 is the reverse complement of the molecule's right-hand continuation, in
reversed order).  We assert the *actual* captured literals and the trimmed
insert, not merely that a rename template happens to reference a capture.

Barcodes / UMIs are fixed string literals (deterministic), so captured values
are asserted exactly.  For SHARE-Seq the barcode/UMI are read off the bottom
strand on R2, so the expected captured values are the ``rc()`` of the top-strand
literals.
"""

import sys
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # Python < 3.11

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cutseq.run import CutadaptConfig, _build_scheme, _scheme_modifiers  # noqa: E402
import cutseq.grammar as grammar  # noqa: E402
from dnaio import SequenceRecord  # noqa: E402
from cutadapt.info import ModificationInfo  # noqa: E402


RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(s):
    return str(s).translate(RC_TABLE)[::-1]


_ADAPTERS = tomllib.loads(
    (Path(__file__).resolve().parent.parent / "cutseq" / "adapters.toml")
    .read_text(encoding="utf-8")
)


# --- helpers ---------------------------------------------------------------


def _run_paired(name, r1_seq, r2_seq, name_format):
    """Run the real paired modifiers for a built-in scheme on one read pair.

    Reads are built to the physical geometry of the fixture; returns the two
    (trimmed, renamed) records.  Uses the scheme straight from adapters.toml so
    the test tracks the source-of-truth correction.
    """
    scheme = _ADAPTERS[name]["scheme"]
    s = CutadaptConfig()
    s.auto_inline = False
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=True, settings=s)
    mods[-1] = cs.renamer(paired=True, name_format=name_format)
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        r1 = SequenceRecord("x/1", r1_seq, "I" * len(r1_seq))
        r2 = SequenceRecord("x/2", r2_seq, "I" * len(r2_seq))
        i1, i2 = ModificationInfo(r1), ModificationInfo(r2)
        for step in mods:
            if isinstance(step, tuple):
                m1, m2 = step
                n1 = m1(r1, i1) if m1 else r1
                n2 = m2(r2, i2) if m2 else r2
                r1, r2 = n1 or r1, n2 or r2
            else:
                step(r1, r2, i1, i2)
        return r1, r2
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


# --- fixed deterministic literals ------------------------------------------

# Smart-seq3
_SS3_R1SITE = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
_SS3_TAG = "ATTGCGCAATG"
_SS3_UMI = "ACGTACGT"            # 8-nt UMI
_SS3_R2SITE = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
_INS30 = "GATTACAGACTTACAGACTTACAGACTTAC"

# 10x GEX half (used by 10X_RNA_ATAC / 10X_RNA_ATAC_MULTI).  The upstream
# TruSeq R1 site is TCTTTCCCTACACGACGCTCTTCCGATCT and anneals upstream of the
# read, so R1 starts after it (at the cell barcode).
_10X_BC = "ACGTACGTACGTACGT"      # 16-nt cell barcode
_10X_UMI = "ACGTACGTACGT"         # 12-nt UMI
_10X_INS = "CGATCGTACGATCGTACGATCGTACGA"

# SHARE-Seq (RNA + ATAC)
_MOME = "AGATGTGTATAAGAGACAG"     # s5 ME (i5-side)
_ME_R2 = "CTGTCTCTTATACACATCT"    # ME R2 (i7-side)
_S7 = "CCGAGCCCACGAGAC"
_R1L = "TCGGACGATCATGGG"
_R1P = "CAAGTATGCAGCGCG"
_R2L = "CTCAAGCACGTGGAT"
_R2P = "AGTCGTACGCCGATG"
_R3L = "CGAAACATCGGCCAC"
_SHARE_UMI = "TGATAACCGA"          # 10-nt mRNA-only UMI
_SHARE_B1 = "TACTCGAC"
_SHARE_B2 = "ATCCGTCA"
_SHARE_B3 = "CGACCGGC"
_SHARE_INS = "TAGGCGTCGATGCCGATCCCACGGA"


# --- Smart-seq3 ------------------------------------------------------------


def test_smartseq3_tag_umi_ggg_insert():
    """Smart-seq3 R1 reads: R1 site + t7-tag + 8-nt UMI + poly-G + insert.

    The R2 site (TruSeq R2) is a real read-2 sequencing SITE, so it anneals
    upstream and is NOT part of R2; a short insert makes R2 = rc(insert).  The
    8-nt UMI is on the left (R1) arm and is captured as-is (top strand).
    """
    r1 = _SS3_R1SITE + _SS3_TAG + _SS3_UMI + "GGG" + _INS30
    r2 = _rc(_INS30)
    out1, out2 = _run_paired("SMART_SEQ3", r1, r2, "{id}_umi:{1}")
    assert out1.sequence == _INS30                     # insert preserved
    assert out2.sequence == _rc(_INS30)               # insert preserved (rc)
    assert out1.name == f"x/1_umi:{_SS3_UMI}"         # UMI captured from R1


def test_smartseq3_internal_no_tag_umi_ggg():
    """Smart-seq3 internal fragment: insert directly flanked by the R1 and R2
    sites, with no t7-tag / UMI / poly-G.  No captures; insert preserved."""
    r1 = _SS3_R1SITE + _INS30
    r2 = _rc(_INS30)
    out1, out2 = _run_paired("SMART_SEQ3_INTERNAL", r1, r2, "{id}")
    assert out1.sequence == _INS30
    assert out2.sequence == _rc(_INS30)
    assert out1.name == "x/1"
    assert out2.name == "x/2"


# --- 10x Multiome GEX half -------------------------------------------------


def test_10x_rna_atac_gex_half_cell_bc_umi_insert():
    """10X_RNA_ATAC is corrected to the GEX (RNA) half only: the ATAC and GEX
    libraries are two separate libraries, not one concatenated molecule.  This
    GEX half equals 10X_RNA_V3: 16-nt cell barcode + 12-nt UMI + insert, with
    the TruSeq R1 site annealed UPSTREAM (so the R1 read starts after it).
    """
    # R1 begins after the upstream TruSeq R1 site: cell barcode + UMI + insert.
    r1 = _10X_BC + _10X_UMI + _10X_INS
    r2 = _rc(_10X_INS)
    out1, out2 = _run_paired(
        "10X_RNA_ATAC", r1, r2, "{id}_bc:{1}_umi:{2}"
    )
    assert out1.sequence == _10X_INS
    assert out2.sequence == _rc(_10X_INS)
    assert out1.name == f"x/1_bc:{_10X_BC}_umi:{_10X_UMI}"


def test_10x_rna_atac_multi_gex_half_cell_bc_umi_insert():
    """10X_RNA_ATAC_MULTI is corrected to the GEX half only; the ATAC and
    MULTI-seq libraries are separate libraries.  Behavior equals 10X_RNA_V3."""
    r1 = _10X_BC + _10X_UMI + _10X_INS
    r2 = _rc(_10X_INS)
    out1, out2 = _run_paired(
        "10X_RNA_ATAC_MULTI", r1, r2, "{id}_bc:{1}_umi:{2}"
    )
    assert out1.sequence == _10X_INS
    assert out2.sequence == _rc(_10X_INS)
    assert out1.name == f"x/1_bc:{_10X_BC}_umi:{_10X_UMI}"


# --- SHARE-Seq -------------------------------------------------------------


def _share_r2(umi=None):
    """The physically-correct R2 read for the SHARE-Seq cassette.

    Top strand 5'->3' (after the s5 ME + insert): [mRNA-only UMI] | r1link | BC1
    | (r1' + r2link) | BC2 | (r2' + r3link) | BC3 | ME R2 + s7.  The ME R2 + s7
    is the outermost read-2 SITE, so it anneals upstream and is NOT part of the
    read.  R2 reads the bottom strand from that site inward: rc(BC3) rc(r2'+r3link)
    rc(BC2) rc(r1'+r2link) rc(BC1) rc(r1link) [rc(UMI)] rc(insert) rc(MOME).
    """
    head = (
        _rc(_SHARE_B3) + _rc(_R2P + _R3L) + _rc(_SHARE_B2) + _rc(_R1P + _R2L)
        + _rc(_SHARE_B1) + _rc(_R1L)
    )
    if umi is not None:
        head += _rc(umi)
    return head + _rc(_SHARE_INS) + _rc(_MOME)


def test_share_seq_rna_cassette_umi_and_barcodes_from_r2():
    """SHARE-Seq RNA cassette (mRNA-only 10-nt UMI).  R1 = s5 ME + insert; the
    UMI and the three 8-nt cell barcodes live on the right arm and are read off
    the bottom strand on R2, so their captured values are ``rc()`` of the top-
    strand literals.  Insert preserved on both reads."""
    r1 = _MOME + _SHARE_INS
    r2 = _share_r2(umi=_SHARE_UMI)
    out1, out2 = _run_paired(
        "SHARE_SEQ", r1, r2, "{id}_umi:{1}_b1:{2}_b2:{3}_b3:{4}"
    )
    assert out1.sequence == _SHARE_INS
    assert out2.sequence == _rc(_SHARE_INS)
    assert out1.name == (
        f"x/1_umi:{_rc(_SHARE_UMI)}_b1:{_rc(_SHARE_B1)}_b2:{_rc(_SHARE_B2)}"
        f"_b3:{_rc(_SHARE_B3)}"
    )


def test_share_seq_atac_cassette_no_umi():
    """SHARE-Seq ATAC cassette: identical to the RNA cassette but with NO UMI.
    The three cell barcodes are captured on R2 as ``rc()`` of the top-strand
    literals; insert preserved on both reads."""
    r1 = _MOME + _SHARE_INS
    r2 = _share_r2(umi=None)
    out1, out2 = _run_paired(
        "SHARE_SEQ_ATAC", r1, r2, "{id}_b1:{1}_b2:{2}_b3:{3}"
    )
    assert out1.sequence == _SHARE_INS
    assert out2.sequence == _rc(_SHARE_INS)
    assert out1.name == (
        f"x/1_b1:{_rc(_SHARE_B1)}_b2:{_rc(_SHARE_B2)}_b3:{_rc(_SHARE_B3)}"
    )


# --- write-domain reference sanity ------------------------------------------


@pytest.mark.parametrize(
    "name,rename,n_caps",
    [
        ("SMART_SEQ3", "{id}_umi:{1}", 1),
        ("SMART_SEQ3_INTERNAL", "{id}", 0),
        ("10X_RNA_ATAC", "{id}_rna-cell_barcode:{1}_rna-umi:{2}", 2),
        ("10X_RNA_ATAC_MULTI", "{id}_rna-cell_barcode:{1}_rna-umi:{2}", 2),
        ("SHARE_SEQ", "{id}_umi:{1}_cell-barcode-1:{2}_cell-barcode-2:{3}"
        "_cell-barcode-3:{4}", 4),
        ("SHARE_SEQ_ATAC", "{id}_cell-barcode-1:{1}_cell-barcode-2:{2}"
        "_cell-barcode-3:{3}", 3),
    ],
)
def test_corrected_scheme_recommended_rename_resolves(name, rename, n_caps):
    """The corrected schemes' recommended_rename references only exist in the
    capture set (write-domain sanity, not a runtime claim)."""
    import re

    from cutseq.grammar import parse_scheme

    assert _ADAPTERS[name]["recommended_rename"] == rename
    _, left, right = parse_scheme(_ADAPTERS[name]["scheme"], auto_inline=False)
    ncaps = len([t for t in left + right if t.kind in ("capture", "inline")])
    assert ncaps == n_caps
    for ref in re.findall(r"\{(\d+)\}", rename):
        assert 1 <= int(ref) <= ncaps
