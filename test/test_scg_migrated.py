#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Migrated SCG/seqspec library schemes (adapters.toml).

Every scheme auto-imported from ``scg_to_cutseq.py`` is fully typed (adapter /
barcode / UMI roles are explicit in the grammar) and therefore registered with
``auto_inline = false`` so cutseq's inline-barcode heuristic never re-labels a
structural linker. These tests split into two honest kinds:

* **Static scheme coverage** — every migrated scheme parses, the write-domain
  ``recommended_rename`` resolves only to captures that exist, and the nine
  maps that today contain >=2 apparent sequencing sites on one arm (so a single
  read has no single 5' read-start) fail actionably as ``ambiguous`` while every
  other map compiles.  Whether those nine are genuine multi-site single
  molecules or importer-*flattened* *alternative* modalities is flagged for the
  seqspec catalog correction — we only pin the current observable behavior.

* **Physically-annotated behavioral fixtures** — a handful of curated, well-
  understood methods are exercised through the real modifier pipeline with a
  hand-built molecular read pair whose geometry we know independently (the
  sequencing primers sit UPSTREAM and are absent from the reads). We assert the
  actually-captured barcode/UMI literals and the trimmed insert, not that the
  rename template merely references a capture.

We deliberately do NOT call the private ``readplan._resolve_boundary`` to build
any expectation here — that would re-derive the source's own decision with the
source's own resolver (circular).  Expected boundary facts come from the scheme
molecular sequence and the curated primer DB, or from the constructed fixture.
"""

import random
import re
import sys
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # Python < 3.11

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cutseq.common import load_adapters_no_auto_inline  # noqa: E402
from cutseq.grammar import parse_scheme  # noqa: E402
from cutseq.run import CutadaptConfig, _build_scheme, _scheme_modifiers  # noqa: E402

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


def migrated_names():
    return {m for m, cfg in _ADAPTERS.items()
            if isinstance(cfg, dict) and "scheme" in cfg
            and cfg.get("auto_inline") is False}


def _capture_count(scheme):
    """Number of capture/inline parts declared by a scheme (written order)."""
    _, left, right = parse_scheme(scheme, auto_inline=False)
    return len([t for t in left + right if t.kind in ("capture", "inline")])


def _compile_paired(scheme):
    """Compile a scheme for a paired run; returns the modifier list (raises on
    an unresolvable / ambiguous read-plan map)."""
    s = CutadaptConfig()
    s.auto_inline = False
    cs = _build_scheme(scheme, s)
    return _scheme_modifiers(cs, paired=True, settings=s)


def _run_paired(scheme, r1_seq, r2_seq, name_format):
    """Run the compiled paired modifiers (real pipeline steps) on one synthetic
    read pair.  Reads are built to the physical geometry of the fixture; the
    output name and both trimmed sequences are returned and asserted on.

    Returns ``(out_r1, out_r2)``.
    """
    import cutseq.grammar as grammar

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


def _made(seed):
    return random.Random(seed)


# --- static scheme coverage ------------------------------------------------


def test_all_migrated_schemes_parse_without_auto_inline():
    names = migrated_names()
    assert names, "no migrated schemes flagged auto_inline=false"
    for name in names:
        parse_scheme(_ADAPTERS[name]["scheme"], auto_inline=False)


def test_migrated_names_report_no_auto_inline():
    migrated = migrated_names()
    reported = set(load_adapters_no_auto_inline())
    assert migrated <= reported


def test_migrated_recommended_rename_references_existing_captures():
    """Every ``{N}`` in a ``recommended_rename`` resolves to an existing
    capture/inline part (1-based, written order).  This checks the WRITE-domain
    reference is consistent with the scheme — it does NOT claim the actual
    runtime capture value is correct."""
    for name in sorted(migrated_names()):
        cfg = _ADAPTERS[name]
        rename = cfg.get("recommended_rename")
        if not rename:
            continue
        ncap = _capture_count(cfg["scheme"])
        for ref in re.findall(r"\{(\d+)\}", rename):
            n = int(ref)
            assert 1 <= n <= ncap, (
                f"{name}: recommended_rename references {{{n}}} but the scheme "
                f"declares only {ncap} capture/inline part(s)"
            )


# --- current ambiguous maps (PENDING catalog correction) ------------------
#
# The remaining migrated maps still fail to resolve a single read-start
# boundary.  Each concatenates >=2 apparent sequencing-site records onto one
# arm, so the read-plan regenerates two different starts and fails actionably
# rather than pick one arbitrarily.  The sites below are read straight from the
# scheme's molecular sequence (NOT by querying ``_resolve_boundary`` at test
# time).
#
# Note: 10X_RNA_ATAC, 10X_RNA_ATAC_MULTI and SMART_SEQ3 were previously in the
# ambiguous set because the importer had *flattened* alternative modalities
# (10x Multiome GEX + ATAC; SMART-Seq3's R1/R2/R1 concatenation) into a single
# map.  The seqspec correction re-split them into single-library schemes that
# now compile unambiguously (verified in test_migrated_unambiguous_maps_compile),
# so they were removed here.
#
CURRENTLY_AMBIGUOUS = {
    # 10x Flex: 3' (RNA + surface-protein).  Right arm repeats the 10x R2 site
    # around the protein reads (outer + after the protein UMI linker).
    "10XFB_3PRIME": ("10x-R2(outer)", "10x-R2(inner)", "Nextera-R1"),
    # 10x Flex: 5' (RNA + VDJ).  Right arm has a 10x-R2, a TruSeq-R1 and a
    # Nextera-R2 site.
    "10XFB_5PRIME": ("10x-R2", "TruSeq-R1", "Nextera-R2"),
    "10XFB_VDJ_5PRIME": ("10x-R2", "TruSeq-R1", "Nextera-R2"),
    # ISSAAC-seq (RNA + ATAC in one map).  Right arm has a 10x-R2 and a
    # Nextera-R2 site.
    "ISSAAC_SEQ": ("10x-R2", "Nextera-R2", "Nextera-R1"),
    # SN-M3C-seq / MCT-seq: the LEFT arm interleaves an R1 site and an R2 site
    # (twice), so it has two distinct R1 (and R2) boundaries on one arm.
    "SN_M3C_SEQ": ("R1-site(a)", "R2-site(a)", "R1-site(b)", "R2-site(b)"),
    "SNMCTSEQ": ("R1-site(a)", "R2-site(a)", "R1-site(b)"),
}


def test_migrated_maps_with_multiple_sites_fail_as_ambiguous():
    """The remaining migrated maps below fail to resolve a single read-start
    boundary and raise the *actionable* ``ambiguous`` error.

    We assert only the observable, honest behavior — the source fails loudly
    (``ValueError ... ambiguous``) rather than silently picking one read-start.
    These maps are genuinely multi-site single molecules (10x Flex, ISSAAC-seq,
    SN-M3C/MCT-seq) that the seqspec catalog has not yet re-split; they must not
    be parked with a blanket ``xfail`` that catches any error, and they must not
    be claimed to compile (they do not).

    The 10x Multiome (GEX+ATAC) and SMART-Seq3 maps were re-split into single-
    library schemes by the correction, so they were moved out of this set and now
    compile unambiguously (asserted in ``test_migrated_unambiguous_maps_compile``).
    """
    # none of these names may drift out of the migrated set unnoticed
    assert set(CURRENTLY_AMBIGUOUS) <= migrated_names()
    for name in sorted(CURRENTLY_AMBIGUOUS):
        with pytest.raises(ValueError) as ei:
            _compile_paired(_ADAPTERS[name]["scheme"])
        assert "ambig" in str(ei.value).lower(), name


def test_migrated_unambiguous_maps_compile():
    """Every migrated map NOT in the documented multi-assay set must compile
    (no read-plan ambiguity) for a paired run."""
    for name in sorted(migrated_names()):
        if name in CURRENTLY_AMBIGUOUS:
            continue
        _compile_paired(_ADAPTERS[name]["scheme"])  # must not raise


# --- curated physically-annotated behavioral fixtures ----------------------
#
# The geometry of these fixtures is known independently and stated next to each
# test.  The sequencing primers anneal UPSTREAM of each read, so they are NOT
# part of the read and never appear in the synthesized R1/R2.


def test_10x_rna_v3_molecule_cell_bc_umi_insert():
    """10x-RNA-v3: R1 site + 16-nt cell barcode + 12-nt UMI + insert ; R2 site.

    Read plan: R1 starts after the R1 site -> ``R1 = BC + UMI + insert``;
    R2 reads the bottom strand from the R2 site inward -> ``R2 = rc(insert)``.
    The pipeline must trim to the insert and capture the BC/UMI literals.
    """
    scheme = _ADAPTERS["10X_RNA_V3"]["scheme"]
    rnd = _made(12345)
    bc = "".join(rnd.choice("ACGT") for _ in range(16))
    umi = "".join(rnd.choice("ACGT") for _ in range(12))
    ins = "".join(rnd.choice("ACGT") for _ in range(30))
    r1, r2 = _run_paired(scheme, bc + umi + ins, _rc(ins),
                         "{id}_bc:{1}_umi:{2}")
    assert r1.sequence == ins
    assert r2.sequence == _rc(ins)
    assert r1.name == f"x/1_bc:{bc}_umi:{umi}"


def test_10x_rna_v2_molecule_cell_bc_umi_insert():
    """10x-RNA-v2: 16-nt cell barcode + 10-nt UMI (+ R1 site upstream)."""
    scheme = _ADAPTERS["10X_RNA_V2"]["scheme"]
    rnd = _made(12345)
    bc = "".join(rnd.choice("ACGT") for _ in range(16))
    umi = "".join(rnd.choice("ACGT") for _ in range(10))
    ins = "".join(rnd.choice("ACGT") for _ in range(30))
    r1, r2 = _run_paired(scheme, bc + umi + ins, _rc(ins),
                         "{id}_bc:{1}_umi:{2}")
    assert r1.sequence == ins
    assert r2.sequence == _rc(ins)
    assert r1.name == f"x/1_bc:{bc}_umi:{umi}"


def test_smartseq2_insert_only_reads_untouched():
    """SMART-Seq2 is primer:primer (Nextera-R1 : Nextera-R2).  Both primers sit
    upstream, so with a short read the reads are insert-only: ``R1 = insert``,
    ``R2 = rc(insert)``, untouched by the pipeline (no scaffold to trim)."""
    scheme = _ADAPTERS["SMART_SEQ2"]["scheme"]
    rnd = _made(77)
    ins = "".join(rnd.choice("ACGT") for _ in range(40))
    r1, r2 = _run_paired(scheme, ins, _rc(ins), "{id}")
    assert r1.sequence == ins
    assert r2.sequence == _rc(ins)
    assert r1.name == "x/1"
    assert r2.name == "x/2"


def test_smartseq2_full_opposite_readthrough_trims_adapters():
    """SMART-Seq2 full read-through: R1 reads through to the Nextera-R2 site
    (``insert + Nextera-R2``) and R2 reads through to the Nextera-R1 site
    (``rc(insert) + rc(Nextera-R1)``).  The read-through gates must trim the
    opposite-site adapters and preserve the insert."""
    scheme = _ADAPTERS["SMART_SEQ2"]["scheme"]
    nextera_r1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    nextera_r2 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    rnd = _made(99)
    ins = "".join(rnd.choice("ACGT") for _ in range(40))
    r1 = ins + nextera_r2
    r2 = _rc(ins) + _rc(nextera_r1)
    out1, out2 = _run_paired(scheme, r1, r2, "{id}")
    assert out1.sequence == ins
    assert out2.sequence == _rc(ins)


# --- static read-2 primer biology ------------------------------------------
#
# For each canonical paired-end method, the scheme's OUTERMOST right adapter is
# the 3' sequencing-SITE; the actual read-2 sequencing primer (the oligo the
# instrument anneals) is ``rc(site)``.  Verified straight from the scheme's
# molecular sequence against the curated primer DB — never by calling the
# private boundary resolver.

_R2_READ2_PRIMER = {
    # (scheme outermost-right adapter) -> rc == this read-2 primer
    "10X_RNA_V3": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "10X_RNA_V2": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "DROP_SEQ": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "SMART_SEQ2": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "CEL_SEQ2": "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
}


def test_migrated_read2_primer_is_rc_of_right_outer_adapter():
    """The read-2 sequencing primer is the reverse complement of the scheme's
    outermost right adapter (the 3' sequencing-site).  This is an independent
    molecular fact about the scheme, not re-derived from the resolver."""
    for name, read2 in _R2_READ2_PRIMER.items():
        scheme = _ADAPTERS[name]["scheme"]
        _, _, right = parse_scheme(scheme, auto_inline=False)
        adps = [t.value for t in right if t.kind == "adp"]
        assert adps, name
        outer = adps[-1]  # outermost 3' adapter (closest to p7 / read-2 side)
        assert _rc(outer) == read2, name
