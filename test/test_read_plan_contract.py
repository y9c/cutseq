"""Behavioral regression tests for the *primer-aware read-plan* cutseq fix.

These tests encode the CONTRACT for a read plan that is derived from the
physical molecule, not from any particular compiled-step implementation. The
synthetic scheme under test is::

    ACACGACGCTCTTCCGATCTXXT...T-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC

Interpreted as a top-strand molecular map, written 5' -> 3'::

    [p5 binder ACACGACGCTCTTCCGATCT][XX2][T-run]  <cDNA insert>  [mask8][topUMI8][p7 AGATCGGAAGAGCACACGTC]

Physical layout (what each mate actually reads):

  * R1 reads the top strand 5' -> 3'.  The p5 sequencing *binding site*
    (``ACACGACGCTCTTCCGATCT``) is upstream of the read and is ABSENT from the
    R1 FASTQ start, so R1 starts at the 5' mask ``XX2``::

        R1 raw  = XX2  T-run  insert  [mask8 topUMI8 p7   if full readthrough]

  * R2 is the reverse complement of the top strand.  The p7 binding site
    (``AGATCGGAAGAGCACACGTC``) is upstream of R2 and is ABSENT from the R2
    FASTQ start, so R2 starts at the UMI::

        R2 raw  = rc(topUMI8)  rc(mask8)  rc(insert)
                  [A-run rc(XX2) rc(p5-primer)   if full readthrough]

Contract the fix must satisfy (derived from the molecule, not the current code):

  * R1's 5' mask ``XX2`` and the leading poly-T run are trimmed UNCONDITIONALLY
    (they are at the defined read start, not gated on any preceding adapter).
  * R2's 5' ``rc(topUMI8)`` (captured UMI) and ``rc(mask8)`` are trimmed
    UNCONDITIONALLY.
  * The mirrored 3' read-through cuts (``mask8/topUMI8/p7`` on R1,
    ``A-run/rc(XX2)/rc(p5)`` on R2) fire ONLY when the opposite arm really was
    sequenced (adapter match), never because a read happens to be > 50 bp.
  * Read-through trims never steal bases from the insert, even for long reads.
  * The leading poly-T mirrors to a trailing A-run on R2, and the poly trim is
    END-ANCHORED: interior homopolymer runs inside the insert are preserved.
  * ``{1}`` (the only capture, on the R2/right side) resolves to the UMI AS
    READ on R2 = ``rc(topUMI8)``, identical on both mates.
  * ``{r1.1}`` is the top-strand UMI (``topUMI8``), only when R1 read through;
    ``{r2.1}`` is the R2 raw UMI (``rc(topUMI8)``), whenever complete.
  * No automatic canonicalization is required; explicit ``rc()`` works.

The UMI (``ATGCTACG``) is deliberately non-palindromic (``rc = CGTAGCAT``), so
a strand mix-up on either mate is caught immediately.

These tests encode the CONTRACT for the fix and are written against the
physical molecule, NOT against the current compiled output. Several of them
FAIL on the current (pre-fix) source by design — that is the point; they are the
baseline the source writer must satisfy. No test is written to match the current
implementation. The parent has since defined the incomplete-capture policy
(``{1}``/``{r2.1}`` are empty when a capture is not complete — see the
``test_incomplete_r2_umi_*`` case), so it IS asserted below, not left open.

A section of *adversarial* cases follows (see the appended tests) that encode
the same contract at the edges the owner flagged in the rejected first source:
embedded custom primers inside a merged token (no whole-token skip), a
recognised primer inside a token with an upstream offset, single-end starting
downstream, an upstream ``X3`` not re-cut on the opposite read-through, custom
primers that are ambiguous or absent (actionable failure), an incomplete R2 UMI
not advertised, a failed in-read scaffold not letting the UMI eat the insert, a
tiny insert, and YAML primer + inline + capture labels.  These may FAIL on the
current partly-implemented / mid-edit source; that is the point — they expose
the remaining gaps.
"""

import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CUTSEQ = str(ROOT / ".venv" / "bin" / "cutseq")

sys.path.insert(0, str(ROOT))
from cutseq import grammar  # noqa: E402
from cutseq.run import CutadaptConfig, _build_scheme, _scheme_modifiers  # noqa: E402
from cutseq.primers import is_known_primer  # noqa: E402
from dnaio import SequenceRecord  # noqa: E402
from cutadapt.info import ModificationInfo  # noqa: E402


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# --- contract scheme + physical molecule ------------------------------------

SCHEME = "ACACGACGCTCTTCCGATCTXXT...T-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"

PRIMER1 = "ACACGACGCTCTTCCGATCT"  # p5-side sequencing primer / binding site (absent from R1)
XX2 = "GG"                        # 5' mask, 2 nt
TRUN = "TTTTT"                    # left poly-T run (mirrors to A-run on R2)
MASK8 = "GCATCGTA"                # 3' mask, 8 nt
UMI8 = "ATGCTACG"                 # top-strand UMI (non-palindromic: rc = CGTAGCAT)
P7 = "AGATCGGAAGAGCACACGTC"       # p7 sequencing primer / binding site (absent from R2)

R2_UMI = _rc(UMI8)                # CGTAGCAT  (the UMI as read on R2)

# Deterministic, diverse inserts that avoid adapter / primer fragments and
# long homopolymer runs, so short random adapter hits cannot fire.
INS = "GATCCAGTTTAACGCCTTGGATTGAGAGACAA"
# Contains an interior T-run (preserved on R1) and its rc has an interior A-run
# (preserved on R2) — used to prove the poly trim is end-anchored.
INS_PT = "CTGCGTAGCTCTCAGTGTTGTCCAACCTTTATTGAGCT"
# Long insert (> 50 bp) with the same clean composition, for the "never steals"
# guard on reads longer than ``--force-trim-min-length``.
INS_L = (
    "GATCCAGTTTAACGCCTTGGATTGAGAGACAA"
    "CTGCGTAGCTCTCAGTGTTGTCCAACCTTTATTGAGCT"
    "AGCTCAATAAAGGTTGGACAACACTGAGAGCTACGCAG"
    "GATCCAGTTTAACGCCTTGGATTGAGAGACAA"
)

RC_INS = _rc(INS)
RC_INS_PT = _rc(INS_PT)
RC_INS_L = _rc(INS_L)

# RC of the p5-side primer (the read-through tail on R2).
RC_P5 = _rc(PRIMER1)
# The R2 sequencing primer oligo is the actual read-2 primer, 5' -> 3': for
# Illumina it is the reverse complement of the top-strand right binding site,
# so ``--r2-primer`` takes rc(P7), not P7.
R2_PRIMER_OLIGO = _rc(P7)

# --- variant schemes ---------------------------------------------------------

# Same molecule WITHOUT the leading poly-T run, to isolate the primer-boundary
# (XX2 / UMI+mask unconditional cuts) behaviour from the poly-run handling.
NOPOLY_SCHEME = "ACACGACGCTCTTCCGATCTXX-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"
# An ``X3`` upstream of the left known primer: it lies outside the read (further
# 5' than the primer) and must be EXCLUDED, not consumed, when the read start is
# inferred at the primer.
X3_SCHEME = "XXXACACGACGCTCTTCCGATCTXX-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"

# Actual CUSTOM sequencing primers that are NOT in the built-in primer DB.  The
# R1 primer (P1) is the top-strand 5' site; the R2 primer (Q2) is the actual R2
# oligo, whose top-strand molecular site is rc(Q2) on the right.
P1_CUSTOM = "CGTGCGATCGTACGGTACCGGTCAG"
Q2_CUSTOM = "TAGCGTACGTACCGGATCGTACGT"
CUSTOM_SCHEME = P1_CUSTOM + "XXT...T-XXXXXXXXNNNNNNNN" + _rc(Q2_CUSTOM)

# --- adversarial primer-plan cases (approved-design gaps) --------------------

# CASE 1: a custom primer is EMBEDDED inside a single contiguous uppercase
# token, with an upstream offset-prefix (before it) and a downstream literal
# scaffold (>= 10 nt) after it.  With auto_inline OFF the token is not split.
CUSTOM_UP_L = "GGATCCGA"
CUSTOM_SCAFFOLD_L = "GTCAGGATCCGTCAGTCGATCGTAC"        # 24 nt, visible residual
CUSTOM_UP_R = "CACCGATT"
CUSTOM_SCAFFOLD_R = "GTATCCGATGCAGGTACCTAGCA"          # 22 nt, visible residual
CUSTOM_MERGED_SCHEME = (
    (CUSTOM_UP_L + P1_CUSTOM + CUSTOM_SCAFFOLD_L)
    + "-"
    + (CUSTOM_SCAFFOLD_R + _rc(Q2_CUSTOM) + CUSTOM_UP_R)
)

# CASE 3: an ``X3`` upstream of the left known primer (with the poly run), for
# a FULL read-through pair where R2's read-through must not cut that outer X3.
X3_POLY_SCHEME = ("XXX" + PRIMER1 + "XXT...T-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC")

# CASE 2: single-end with a recognised left primer -> read starts downstream.
SE_PRIMER_SCHEME = "ACACGACGCTCTTCCGATCTXXNNNNNN:"
UMI6 = UMI8[:6]  # 6 nt UMI matching the N6 capture

# CASE 4: the SAME custom primer appears twice at two distinct offsets in the
# left arm -> the read-start boundary is ambiguous and must fail actionably.
AMBIGUOUS_R1_SCHEME = ("AAA" + P1_CUSTOM + "GGG" + P1_CUSTOM + "TTT-" + _rc(Q2_CUSTOM))

# CASE 5 / normal fallback: a primer-boundary scheme that does NOT contain the
# custom primer -> the supplied custom primer is absent.
CASE5_SCHEME = "ACACGACGCTCTTCCGATCTXX-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"

# CASE 8 (corrected): an EXPLICIT read-visible literal scaffold (>=10 nt) on the
# R2 own arm BEFORE a downstream X4 mask, so its absence can be detected by
# sequence match (an unknown ``XXXXXXXX`` mask is positional-only).
SCAFFOLD_TOP = "ACGTTCGATCGGTACCGTAACT"          # 22 nt top-strand scaffold
SCAFF_R2 = _rc(SCAFFOLD_TOP)                     # what R2 actually reads (rc)
WRONG_SCAFFOLD = "G" * 22                        # mismatching 22-mer (not SCAFF_R2)
MASK4_READ = "GATC"                              # 4 actual bases for the X4 mask
SCAFFOLD_SCHEME = (PRIMER1 + "XX-" + "XXXX" + SCAFFOLD_TOP + "NNNNNNNN" + P7)


def _r1_nort(insert=INS):
    return XX2 + TRUN + insert
def _r2_nort(insert=INS):
    return R2_UMI + _rc(MASK8) + _rc(insert)
def _r1_full(insert=INS):
    return XX2 + TRUN + insert + MASK8 + UMI8 + P7
def _r2_full(insert=INS):
    return R2_UMI + _rc(MASK8) + _rc(insert) + _rc(TRUN) + _rc(XX2) + _rc(PRIMER1)
# PARTIAL read-through: the opposite arm is cut off mid-adapter (a 10 bp prefix
# of p7 on R1 / rc(p5) on R2) but the adjacent full UMI+mask are present.
def _r1_partial(insert=INS):
    return XX2 + TRUN + insert + MASK8 + UMI8 + P7[:10]
def _r2_partial(insert=INS):
    return R2_UMI + _rc(MASK8) + _rc(insert) + _rc(TRUN) + _rc(XX2) + RC_P5[:10]


# --- manual drive through the compiled pipeline -----------------------------

_HIGH_QUAL_ALPHA = "56789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def _qual(n, seed=0):
    """A distinct, all-high-Phred quality string (>= 20) so the trailing
    ``QualityTrimmer`` never removes bases and we can assert exact slices."""
    return "".join(_HIGH_QUAL_ALPHA[(seed + i) % len(_HIGH_QUAL_ALPHA)] for i in range(n))


def _run_paired(r1_seq, r2_seq, name_format="{id}", q1=None, q2=None, scheme=SCHEME,
                settings=None):
    """Run the compiled paired modifiers on one synthetic pair and return the
    two (possibly shortened/renamed) ``SequenceRecord`` objects.

    Drives the same public pipeline the CLI uses (``_build_scheme`` +
    ``_scheme_modifiers``) so the tests stay behavioural: they assert on
    ``sequence`` / ``qualities`` / ``name``, never on internals.
    """
    s = settings or CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=True, settings=s)
    mods[-1] = cs.renamer(paired=True, name_format=name_format)  # sets the rename flag
    try:
        q1 = q1 or "I" * len(r1_seq)
        q2 = q2 or "I" * len(r2_seq)
        r1 = SequenceRecord("x/1", r1_seq, q1)
        r2 = SequenceRecord("x/2", r2_seq, q2)
        i1, i2 = ModificationInfo(r1), ModificationInfo(r2)
        # Explicit ``is not None`` checks: a modifier may return an EMPTY read
        # (e.g. UMI+mask cut more than the read length) without discarding it, so
        # ``n or r`` would wrongly resurrect the old read (an empty SequenceRecord
        # is falsy). We also honour a paired step's returned pair when present.
        for step in mods:
            if isinstance(step, tuple):
                m1, m2 = step
                n1 = m1(r1, i1) if m1 is not None else r1
                n2 = m2(r2, i2) if m2 is not None else r2
                if n1 is not None:
                    r1 = n1
                if n2 is not None:
                    r2 = n2
            else:
                result = step(r1, r2, i1, i2)
                if result is not None:
                    r1, r2 = result
        return r1, r2
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


def _run_single(read_seq, name_format="{id}", scheme=SE_PRIMER_SCHEME, q=None,
                settings=None):
    """Run the compiled single-end modifiers on one read; return the result."""
    s = settings or CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    mods[-1] = cs.renamer(paired=False, name_format=name_format)  # sets the rename flag
    try:
        read = SequenceRecord("x", read_seq, q or "I" * len(read_seq))
        info = ModificationInfo(read)
        for m in mods:
            out = m(read, info)
            if out is not None:
                read = out
        return read
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


# --- pure no-read-through: trim to the insert -------------------------------

def test_no_readthrough_r1_r2_trim_to_insert():
    """Neither mate reaches the opposite arm: R1 -> insert, R2 -> rc(insert).
    R1's 5' mask is cut unconditionally and R2's UMI/mask are cut
    unconditionally, regardless of any adapter match."""
    r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format="{id}")
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence


def test_no_readthrough_names_resolve_r2_umi():
    """``{1}`` is the single capture (R2/right anchor) -> the UMI as read on
    R2 = rc(topUMI8), identical on both mates even with no read-through."""
    r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format="{id}_{1}")
    assert r1.name == f"x/1_{R2_UMI}", r1.name
    assert r2.name == f"x/2_{R2_UMI}", r2.name


def test_no_readthrough_r1_umi_absent_r2_present():
    """``{r1.1}`` (top-strand UMI) is absent without read-through; ``{r2.1}``
    (R2 raw UMI) is present as soon as R2 reads it."""
    r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format="{id}_a={r1.1}_b={r2.1}")
    assert r1.name == f"x/1_a=_b={R2_UMI}", r1.name
    assert r2.name == f"x/2_a=_b={R2_UMI}", r2.name


# --- full read-through: still trim to the insert ----------------------------

def test_full_readthrough_trims_to_insert():
    """Both mates read all the way through the opposite arm; the read-through
    arms are cut but the insert (R1) / rc(insert) (R2) is preserved."""
    r1, r2 = _run_paired(_r1_full(), _r2_full(), name_format="{id}")
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence


def test_full_readthrough_capture_side_resolution():
    """With read-through, ``{r1.1}`` = top-strand UMI, ``{r2.1}`` = R2 raw UMI,
    and unprefixed ``{1}`` = R2 raw UMI (anchor read)."""
    r1, r2 = _run_paired(
        _r1_full(), _r2_full(),
        name_format="{id}_r1={r1.1}_r2={r2.1}_unp={1}",
    )
    assert r1.name == f"x/1_r1={UMI8}_r2={R2_UMI}_unp={R2_UMI}", r1.name
    assert r2.name == f"x/2_r1={UMI8}_r2={R2_UMI}_unp={R2_UMI}", r2.name


def test_explicit_rc_reaches_mirror_umi():
    """``rc({r2.1})`` is the top-strand UMI and equals ``{r1.1}`` in full
    read-through; no automatic canonicalization is required."""
    r1, _ = _run_paired(
        _r1_full(), _r2_full(),
        name_format="{id}_rc=rc({r2.1})_top={r1.1}",
    )
    assert r1.name == f"x/1_rc={UMI8}_top={UMI8}", r1.name


def test_partial_readthrough_trims_readthrough_arm():
    """The opposite arm is read only as a short prefix (10 bp of p7 on R1, 10 bp
    of rc(p5) on R2) after the complete adjacent UMI+mask; the partial
    read-through is still recognized and trimmed, so each mate -> insert."""
    r1, r2 = _run_paired(_r1_partial(), _r2_partial(), name_format="{id}")
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence


# --- mirror cuts must never steal bases (read > 50) -------------------------

def test_long_read_never_steals_mirror_bases():
    """A 150+ bp read that does NOT contain the opposite arm must not have the
    mirrored read-through cuts steal insert bases just because it is long."""
    r1, r2 = _run_paired(_r1_nort(INS_L), _r2_nort(INS_L), name_format="{id}")
    assert len(r1.sequence) > 50
    assert len(r2.sequence) > 50
    assert r1.sequence == INS_L, r1.sequence
    assert r2.sequence == RC_INS_L, r2.sequence


# --- asymmetric mates -------------------------------------------------------

def test_asymmetric_r1_readthrough_r2_not():
    """R1 reads through the p7 arm, R2 stops before its read-through: each mate
    still trims to its own insert/rc(insert)."""
    r1, r2 = _run_paired(
        _r1_full(), _r2_nort(), name_format="{id}_a={r1.1}_b={r2.1}",
    )
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence
    assert r1.name == f"x/1_a={UMI8}_b={R2_UMI}", r1.name
    assert r2.name == f"x/2_a={UMI8}_b={R2_UMI}", r2.name


def test_asymmetric_r2_readthrough_r1_not():
    """R2 reads through the p5 arm, R1 stops before its read-through: R1's
    top-strand UMI capture stays empty, R2's stays populated."""
    r1, r2 = _run_paired(
        _r1_nort(), _r2_full(), name_format="{id}_a={r1.1}_b={r2.1}_c={1}",
    )
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence
    assert r1.name == f"x/1_a=_b={R2_UMI}_c={R2_UMI}", r1.name
    assert r2.name == f"x/2_a=_b={R2_UMI}_c={R2_UMI}", r2.name


# --- qualities slice exactly ------------------------------------------------

def test_qualities_slice_exactly_with_trim():
    """Permitting the read to survive, each trim slices the quality string with
    the same offsets as the sequence."""
    q1 = _qual(len(_r1_nort()), seed=3)   # covers XX2 + T-run + insert
    q2 = _qual(len(_r2_nort()), seed=7)   # covers UMI + mask + rc(insert)
    r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format="{id}", q1=q1, q2=q2)
    assert r1.sequence == INS
    assert r1.qualities == q1[len(XX2) + len(TRUN):], r1.qualities
    assert r2.sequence == RC_INS
    assert r2.qualities == q2[len(R2_UMI) + len(_rc(MASK8)):], r2.qualities


# --- poly tails end-anchored; interior runs preserved -----------------------

def test_polyt_end_anchored_preserves_interior_t_on_r1():
    """The leading poly-T is trimmed from R1's 5'; an interior T-run inside the
    insert is preserved."""
    r1, r2 = _run_paired(_r1_nort(INS_PT), _r2_nort(INS_PT), name_format="{id}")
    assert r1.sequence == INS_PT, r1.sequence
    assert r2.sequence == RC_INS_PT, r2.sequence
    # sanity: the insert really does contain an interior T-run (and rc has an A-run)
    assert "TT" in INS_PT and "AA" in RC_INS_PT


def test_polyt_mirrors_to_a_run_preserves_interior_a_on_r2():
    """The leading poly-T mirrors to a trailing A-run on R2 read-through, and it
    is END-anchored so an interior A-run in rc(insert) is not truncated."""
    r1, r2 = _run_paired(_r1_full(INS_PT), _r2_full(INS_PT), name_format="{id}")
    assert r1.sequence == INS_PT, r1.sequence
    assert r2.sequence == RC_INS_PT, r2.sequence
    # rc(INS_PT) contains an interior A-run that must survive the A-run trim.
    assert "AA" in RC_INS_PT


# --- primer-boundary isolation: no-poly variant + upstream-excluded token ---

def test_nopoly_variant_no_readthrough_isolates_primer_boundary():
    """Without the poly-T run the contract still holds: XX2 (R1) and UMI/mask
    (R2) are cut unconditionally, so no-read-through -> insert / rc(insert)."""
    r1, r2 = _run_paired(XX2 + INS, R2_UMI + _rc(MASK8) + RC_INS,
                         name_format="{id}", scheme=NOPOLY_SCHEME)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence


def test_nopoly_variant_full_readthrough_isolates_primer_boundary():
    """No-poly full read-through: the whole opposite arm is trimmed to the
    insert, showing the 5' cuts and the read-through work independently of the
    poly-run handling."""
    r1 = XX2 + INS + MASK8 + UMI8 + P7
    r2 = R2_UMI + _rc(MASK8) + RC_INS + _rc(XX2) + RC_P5
    r1o, r2o = _run_paired(r1, r2, name_format="{id}", scheme=NOPOLY_SCHEME)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence


def test_upstream_excluded_token_before_known_primer():
    """An ``X3`` placed upstream (5') of the left known primer lies outside the
    read. When the read start is inferred at the primer, the outer ``X3`` must be
    EXCLUDED, not consumed (its bases must not leak into the trimmed insert)."""
    r1, r2 = _run_paired(XX2 + INS, R2_UMI + _rc(MASK8) + RC_INS,
                         name_format="{id}", scheme=X3_SCHEME)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence


# --- custom sequencing primers located in the scheme --------------------------

def test_custom_primers_in_scheme_no_readthrough_via_cli():
    """Actual custom primers (not in the DB) appear in the scheme as outer
    adapters; ``--r1-primer``/``--r2-primer`` give their true 5'->3' oligo
    orientation (R1 = P1, R2 = Q2, right molecular site = rc(Q2)) and locate the
    read start, so no-read-through trims to insert."""
    r1 = XX2 + TRUN + INS
    r2 = R2_UMI + _rc(MASK8) + RC_INS
    (n1, s1), (n2, s2) = _cli(
        CUSTOM_SCHEME, r1, r2, "{id}",
        "--r1-primer", P1_CUSTOM, "--r2-primer", Q2_CUSTOM,
    )
    assert s1 == INS and s2 == RC_INS, (s1, s2)
    assert n1 == "@x/1" and n2 == "@x/2", (n1, n2)


def test_custom_primers_in_scheme_full_readthrough_via_cli():
    """Same custom-primer scheme, but both mates read through the opposite arm;
    the read-through arm is trimmed and the insert preserved."""
    r1 = XX2 + TRUN + INS + MASK8 + UMI8 + _rc(Q2_CUSTOM)
    r2 = R2_UMI + _rc(MASK8) + RC_INS + _rc(TRUN) + _rc(XX2) + _rc(P1_CUSTOM)
    (n1, s1), (n2, s2) = _cli(
        CUSTOM_SCHEME, r1, r2, "{id}",
        "--r1-primer", P1_CUSTOM, "--r2-primer", Q2_CUSTOM,
    )
    assert s1 == INS and s2 == RC_INS, (s1, s2)
    assert n1 == "@x/1" and n2 == "@x/2", (n1, n2)


# --- orientation markers are equivalent (paired, no auto-rc) ----------------

def test_orientation_markers_equivalent_trimming():
    """``+``, ``-`` and ``:`` split the library identically for trimming (the
    marker only matters for single-end ``--auto-rc``, which is ignored here)."""
    outs = []
    for marker in ("+", ":", "-"):
        scheme = SCHEME.replace("-", marker, 1) if marker != "-" else SCHEME
        assert marker in scheme
        r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format="{id}", scheme=scheme)
        outs.append((r1.sequence, r2.sequence))
    assert all(o == outs[0] for o in outs), outs
    assert outs[0] == (INS, RC_INS)


# --- non-primer custom scaffold must stay trimmed (DBiT precedent) ----------

def test_no_primer_custom_scaffold_is_trimmed_not_treated_as_binder():
    """A custom outer adapter that is NOT a known sequencing primer must still
    be trimmed from the read (it is part of the read, not an upstream binding
    site). This guards against the read-plan fix blanket-marking outer adapters
    as absent binders (DBiT-seq precedent)."""
    handle = "CAAGCGTTGGCTTCTCGCATCT"          # DBiT / m6A-ARTR handle, not a primer
    assert not is_known_primer(handle)
    scheme = handle + ":"
    r1, _ = _run_paired(handle + INS, _rc(INS), name_format="{id}", scheme=scheme)
    assert r1.sequence == INS, r1.sequence


# --- built-in primer merged with a downstream literal scaffold ---------------

def test_builtin_primer_merged_with_scaffold_preserves_residual():
    """A built-in sequencing primer concatenated with a downstream literal
    scaffold (a single uppercase run) is auto-split into primer + inline; the
    read plan must preserve and trim the scaffold *residual* (no whole-token
    skip), so the read still reaches the insert."""
    merged = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTGGGG-XXXXXXXXNNNNNNNNAGATCGGAAGAGCACACGTC"
    scaffold = "GGGG"
    r1, r2 = _run_paired(scaffold + INS, R2_UMI + _rc(MASK8) + _rc(INS),
                         name_format="{id}", scheme=merged)
    assert r1.sequence == INS, r1.sequence        # residual trimmed, not skipped
    assert r2.sequence == RC_INS, r2.sequence


# --- YAML part scheme drives the same contract ------------------------------

def test_yaml_scheme_label_capture_and_trim():
    """The same molecule expressed as an explicit YAML part map (read-visible)
    trims identically and resolves the labeled capture as the R2 UMI."""
    import yaml
    data = {
        "parts": [
            {"type": "adapter", "seq": PRIMER1},
            {"type": "mask", "length": 2},
            {"type": "polytail", "base": "T"},
            {"type": "insert", "value": "-"},
            {"type": "mask", "length": 8},
            {"type": "capture", "length": 8, "label": "umi"},
            {"type": "adapter", "seq": P7},
        ],
    }
    with tempfile.TemporaryDirectory() as td:
        path = f"{td}/scheme.yaml"
        with open(path, "w") as fh:
            yaml.safe_dump(data, fh)
        r1, r2 = _run_paired(_r1_full(), _r2_full(), name_format="{id}_{umi}", scheme=path)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence
    assert r1.name == f"x/1_{R2_UMI}", r1.name
    assert r2.name == f"x/2_{R2_UMI}", r2.name


# --- explicit back-adapter token does not crash -----------------------------

def test_explicit_back_adapter_token_does_not_crash():
    """A scheme containing a ``>`` back-adapter token compiles and runs the
    pipeline without crashing."""
    scheme = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT:>TCGTCGGCAGCGTC"
    s = CutadaptConfig()
    cs = _build_scheme(scheme, s)
    mods = _scheme_modifiers(cs, paired=False, settings=s)
    mods[-1] = cs.renamer(paired=False, name_format="{id}")
    grammar._RENAME_NEEDS_CAPTURES = True
    try:
        read_seq = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" + "TTTTACGT"
        read = SequenceRecord("y", read_seq, "I" * len(read_seq))
        info = ModificationInfo(read)
        out = read
        for m in mods:
            nxt = m(out, info) if m else out
            if nxt is not None:
                out = nxt
        # The back adapter is 'own'-side: at minimum we must produce a read.
        assert isinstance(out, SequenceRecord)
    finally:
        grammar._RENAME_NEEDS_CAPTURES = False
        grammar._capture_registry.clear()


# --- CLI end-to-end (public binary) ------------------------------------------

def _cli(scheme, r1, r2, name_format, *extra, min_len=5):
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            f1.write(f"@x/1\n{r1}\n+\n{'I' * len(r1)}\n")
            f2.write(f"@x/2\n{r2}\n+\n{'I' * len(r2)}\n")
        p = subprocess.run(
            [CUTSEQ, "-A", scheme, "-O", f"{td}/o", "--rename", name_format,
             "-m", str(min_len), *extra, f"{td}/t1.fq", f"{td}/t2.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr

        def read(side):
            with gzip.open(f"{td}/o_trimmed_{side}.fastq.gz", "rt") as fh:
                lines = [ln.rstrip("\n") for ln in fh]
            return lines[0], lines[1]

        return read("R1"), read("R2")


def test_cli_no_readthrough_trim_and_name():
    n1, s1 = _cli(SCHEME, _r1_nort(), _r2_nort(), "{id}_{1}")[0]
    assert s1 == INS, s1
    assert n1 == f"@x/1_{R2_UMI}", n1


def test_cli_full_readthrough_trim():
    (n1, s1), (n2, s2) = _cli(SCHEME, _r1_full(), _r2_full(), "{id}")
    assert s1 == INS and s2 == RC_INS, (s1, s2)
    assert n1 == "@x/1" and n2 == "@x/2", (n1, n2)


def test_cli_r1_r2_primer_oligo_orientation_accepted():
    """``--r1-primer``/``--r2-primer`` take the actual sequencing oligo 5'->3'
    orientation: --r2-primer is the rc of the top-strand p7 binding site (the
    true read-2 primer), not the binding site itself. They identify upstream
    binding sites and do not inject any 5' trim; reads still trim to insert."""
    n1, s1 = _cli(
        SCHEME, _r1_nort(), _r2_nort(), "{id}",
        "--r1-primer", PRIMER1, "--r2-primer", R2_PRIMER_OLIGO,
    )[0]
    assert s1 == INS, s1
    assert n1 == "@x/1", n1


def _umi_for(k):
    """A deterministic, non-palindromic 8-nt top-strand UMI for record k."""
    v = (k * 7 + 3) % 65536
    n, out = v, ""
    for _ in range(8):
        out = "ACGT"[n % 4] + out
        n //= 4
    if _rc(out) == out:
        out = "A" + out[1:]
    return out


def test_cli_multiprocessing_deterministic():
    """``-t N`` (parallel runner) must not change trimmed sequences or resolved
    capture names, and must not bleed captures between records.  It uses a
    VARYING non-palindromic UMI per record so a cross-record / cross-worker
    registry leak would show up in the names.  Compares EVERY record (all four
    lines each) and verifies the exact expected per-record name + trimmed seq."""
    def pair(k):
        u = _umi_for(k)
        if u == _rc(u):
            u = u[:-1] + ("A" if u[-1] != "A" else "C")
        r1 = XX2 + TRUN + INS
        r2 = _rc(u) + _rc(MASK8) + RC_INS
        return r1, r2, _rc(u)   # (R1 seq, R2 seq, expected {1} = R2 UMI)

    records = [pair(k) for k in range(40)]
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            for k, (r1, r2, _) in enumerate(records):
                f1.write(f"@x{k}/1\n{r1}\n+\n{'I' * len(r1)}\n")
                f2.write(f"@x{k}/2\n{r2}\n+\n{'I' * len(r2)}\n")

        def run_at(threads):
            p = subprocess.run(
                [CUTSEQ, "-A", SCHEME, "-O", f"{td}/o{threads}",
                 "--rename", "{id}_{1}", "-m", "5", "-t", str(threads),
                 f"{td}/t1.fq", f"{td}/t2.fq"],
                capture_output=True, text=True,
            )
            assert p.returncode == 0, p.stderr
            out = {}
            for side in ("R1", "R2"):
                with gzip.open(f"{td}/o{threads}_trimmed_{side}.fastq.gz", "rt") as fh:
                    lines = [ln.rstrip("\n") for ln in fh]
                out[side] = [lines[i:i + 4] for i in range(0, len(lines), 4)]
            return out

        mono, multi = run_at(1), run_at(4)
        # (1) identical records, same order, at every worker count
        assert mono == multi, (mono, multi)
        # (2) exact per-record expected name + trimmed sequence (all 40)
        for side, mate in (("R1", 1), ("R2", 2)):
            for k, rec in enumerate(mono[side]):
                umi = records[k][2]
                assert rec[0] == f"@x{k}/{mate}_{umi}", (side, k, rec[0])
                assert rec[1] == (INS if side == "R1" else RC_INS), (side, k, rec[1])
        # (3) count check
        assert len(mono["R1"]) == 40 and len(mono["R2"]) == 40


def test_polyt_2dot_alias_equivalent():
    """``T..T`` (2 dots) is a requested alias for the ``T...T`` poly-tail, and
    must produce the same read plan as the 3-dot form."""
    scheme2 = SCHEME.replace("T...T", "T..T", 1)
    r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format="{id}", scheme=scheme2)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence


# ============================================================================
# Adversarial cases for the APPROVED primer-aware read-plan design.
# Each encodes the parent/owner contract derived from the physical molecule.
# These may FAIL on the current (partly-implemented / mid-edit) source; that is
# the point — they expose the remaining design gaps.  They are written via the
# public pipeline/CLI only; no internal (readplan) classes are inspected.
# ============================================================================


def test_embedded_custom_primer_preserves_downstream_scaffold():
    """A custom (non-DB) primer embedded INSIDE one merged uppercase token, with
    an upstream offset-prefix and a downstream literal scaffold (>= 10 nt),
    auto_inline off: the upstream offset + primer are NOT consumed, and the
    residual downstream scaffold IS matched (anchored) and trimmed, so each mate
    trims to its insert.  This is the "no whole-token skip" case where the
    primer boundary is interior to a token."""
    r1 = CUSTOM_SCAFFOLD_L + INS
    r2 = _rc(CUSTOM_SCAFFOLD_R) + RC_INS
    (n1, s1), (n2, s2) = _cli(
        CUSTOM_MERGED_SCHEME, r1, r2, "{id}",
        "--no-auto-inline", "--r1-primer", P1_CUSTOM, "--r2-primer", Q2_CUSTOM,
    )
    assert s1 == INS, s1
    assert s2 == RC_INS, s2
    assert n1 == "@x/1", n1
    assert n2 == "@x/2", n2


def test_single_end_known_primer_read_start_downstream():
    """Single-end with a recognised left primer must start DOWNSTREAM of it
    (R1 = XX2 + UMI + insert), not behind a legacy gate on the primer: the R1
    mask and capture are cut unconditionally -> insert, capture = UMI."""
    read = _run_single(XX2 + UMI6 + INS, name_format="{id}_{1}")
    assert read.sequence == INS, read.sequence
    assert read.name == f"x_{UMI6}", read.name


def test_x3_before_primer_full_readthrough_r2_must_not_cut_x3():
    """An ``X3`` upstream of the left known primer is OUTSIDE the read.  In a
    FULL read-through pair, R1 must not consume it and R2's read-through walk
    must NOT process the outer ``X3`` as an inward cut (i.e. it must not loop
    the whole arm past the boundary).  Own-arm / mirror cuts never steal."""
    r1 = XX2 + TRUN + INS + MASK8 + UMI8 + P7
    r2 = R2_UMI + _rc(MASK8) + RC_INS + _rc(TRUN) + _rc(XX2) + RC_P5
    r1o, r2o = _run_paired(r1, r2, name_format="{id}", scheme=X3_POLY_SCHEME)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence


def test_custom_primer_appearing_twice_is_ambiguous():
    """If the SAME custom primer resolves to two distinct read-start boundaries
    in the same arm, the plan must FAIL actionably (ValueError / CLI nonzero)
    rather than guess.  Same-boundary aliases are not the target here."""
    s = CutadaptConfig()
    s.auto_inline = False
    s.r1_primer = P1_CUSTOM
    cs = _build_scheme(AMBIGUOUS_R1_SCHEME, s)
    with pytest.raises(ValueError) as excinfo:
        _scheme_modifiers(cs, paired=True, settings=s)
    assert "ambig" in str(excinfo.value).lower()


def test_custom_primer_supplied_but_absent_fails_actionably():
    """A custom primer explicitly supplied but ABSENT from the scheme is a
    configuration error and must fail actionably, not silently no-op."""
    s = CutadaptConfig()
    s.r1_primer = P1_CUSTOM
    cs = _build_scheme(CASE5_SCHEME, s)
    with pytest.raises(ValueError):
        _scheme_modifiers(cs, paired=True, settings=s)


def test_normal_scheme_without_custom_primer_is_unchanged():
    """A normal primer-boundary scheme with NO supplied custom primer and NO
    visible scaffold must be unchanged (builds fine; the known left primer is
    the boundary)."""
    s = CutadaptConfig()
    cs = _build_scheme(CASE5_SCHEME, s)
    assert cs.orientation == "-"
    assert [(t.kind, t.value) for t in cs.left] == [("adp", PRIMER1), ("mask", 2)]


def test_incomplete_r2_umi_not_advertised_and_r2_empty():
    """A 3 bp R2 that starts only a fraction of the 8 nt UMI reads is NOT
    advertised as a complete capture: ``{1}`` / ``{r2.1}`` stay empty, R2 trims
    to empty, and R1 is unchanged by the other mate's shortness."""
    r1, r2o = _run_paired(XX2 + TRUN + INS, R2_UMI[:3], name_format="{id}_{1}")
    assert r1.sequence == INS, r1.sequence
    assert r2o.sequence == "", r2o.sequence
    assert r1.name == "x/1_", r1.name
    assert r2o.name == "x/2_", r2o.name
    r1b, _ = _run_paired(XX2 + TRUN + INS, R2_UMI[:3], name_format="{id}_{r2.1}")
    assert r1b.name == "x/1_", r1b.name


def test_failed_in_read_scaffold_does_not_let_umi_consume_insert():
    """Corrected: uses an EXPLICIT read-visible literal scaffold (>=10 nt) on the
    R2 own arm BEFORE a downstream X4 mask, so its presence can be detected by
    sequence match (an unknown ``XXXXXXXX`` mask is positional-only and cannot be
    "missing").  When the scaffold MATCHES, UMI + scaffold + X4 are trimmed and
    the insert is preserved; when the scaffold MISMATCHES, the failed sequence
    match BLOCKS the downstream X4 cut so the remainder after the UMI is left
    unchanged (no arbitrary insert is consumed).  R1 is unaffected in both."""
    s = CutadaptConfig()
    s.auto_inline = False

    # Matching: scaffold present as SCAFF_R2 -> UMI + scaffold + X4 all trimmed.
    matching = R2_UMI + SCAFF_R2 + MASK4_READ + RC_INS
    r1, r2o = _run_paired(XX2 + INS, matching, name_format="{id}_{1}",
                          scheme=SCAFFOLD_SCHEME, settings=s)
    assert r1.sequence == INS, r1.sequence
    assert r2o.sequence == RC_INS, r2o.sequence
    assert r1.name == f"x/1_{R2_UMI}", r1.name
    assert r2o.name == f"x/2_{R2_UMI}", r2o.name

    # Mismatch: scaffold region reads a different sequence -> sequence match fails
    # -> the downstream X4 mask is NOT cut; the remainder after UMI is preserved.
    mismatch = R2_UMI + WRONG_SCAFFOLD + MASK4_READ + RC_INS
    r1b, r2b = _run_paired(XX2 + INS, mismatch, name_format="{id}",
                           scheme=SCAFFOLD_SCHEME, settings=s)
    assert r1b.sequence == INS, r1b.sequence
    assert r2b.sequence == WRONG_SCAFFOLD + MASK4_READ + RC_INS, r2b.sequence


def test_tiny_insert_no_poly_full_readthrough_no_steal():
    """A 2 nt insert with no poly run, full read-through: own-arm and mirror cuts
    must not steal from the tiny insert or each other's arm; UMI stays correct."""
    small = "GA"
    r1 = XX2 + small + MASK8 + UMI8 + P7
    r2 = R2_UMI + _rc(MASK8) + _rc(small) + _rc(XX2) + RC_P5
    r1o, r2o = _run_paired(r1, r2, name_format="{id}_{1}", scheme=NOPOLY_SCHEME)
    assert r1o.sequence == small, r1o.sequence
    assert r2o.sequence == _rc(small), r2o.sequence
    assert r1o.name == f"x/1_{R2_UMI}", r1o.name
    assert r2o.name == f"x/2_{R2_UMI}", r2o.name


def test_yaml_primer_inline_capture_labels_and_ensure_inline_barcode():
    """A YAML primer plan with primer + inline (label bc) + capture (label umi)
    plus a right capture (label umi3): the inline barcode at the read start is
    captured as ``{bc}`` (it is the visible read-start scaffold, not skipped),
    and ``--ensure-inline-barcode`` runs without error for a matching read."""
    import yaml

    data = {
        "parts": [
            {"type": "adapter", "seq": PRIMER1},
            {"type": "inline", "seq": "ATCACG", "label": "bc"},
            {"type": "capture", "length": 8, "label": "umi"},
            {"type": "insert", "value": "-"},
            {"type": "mask", "length": 8},
            {"type": "capture", "length": 8, "label": "umi3"},
            {"type": "adapter", "seq": P7},
        ],
    }
    with tempfile.TemporaryDirectory() as td:
        path = f"{td}/s.yaml"
        with open(path, "w") as fh:
            yaml.safe_dump(data, fh)
        r1 = "ATCACG" + UMI8 + INS
        r2 = R2_UMI + _rc(MASK8) + RC_INS
        n1, s1 = _cli(
            path, r1, r2,
            "{id}_bc={bc}_umi={umi}_u3={umi3}",
            "--ensure-inline-barcode",
        )[0]
    assert s1 == INS, s1
    assert n1 == f"@x/1_bc=ATCACG_umi={UMI8}_u3={R2_UMI}", n1

# ============================================================================
# Cases flagged by the parent audit (source writer replacing the exact-string
# executor with cutadapt primitives).  All expectations derive from the physical
# molecule; no internal readplan classes are referenced.
# ============================================================================

# CASE A: full read-through adapter followed by EXTRA non-adapter bases.
EXTRA_TAIL = "ACGT"

# CASE B: no front token at all — primer, insert, primer (DSLIGATION-like).
DSLIGATION_SCHEME = PRIMER1 + ":" + P7

# CASE C: the same BUILTIN (known) primer twice at distinct boundaries -> ambiguity.
BUILTIN_AMBIGUOUS_SCHEME = "AAA" + PRIMER1 + "GGG" + PRIMER1 + "TTT-" + _rc(P7)

# CASE F (single-end) / CASE G (mixed): a right-hand NON-primer read-through
# scaffold that is read-visible, so it is NOT a recognised binding site.
RIGHT_SCAFFOLD = "ACGTTCGATCGGTACCGTAACT"
SE_RIGHT_SCHEME = PRIMER1 + "XX:" + RIGHT_SCAFFOLD
MIXED_SCHEME = PRIMER1 + ":" + RIGHT_SCAFFOLD

# CASE E: an embedded BUILTIN R2 binding site (P7) inside one merged uppercase
# token, with a visible residual prefix (before it) and an upstream suffix
# (after it), auto_inline off.
EMBED_R2_RESIDUAL = "GTATCCGATGCAGGTACCTAGCA"
EMBED_R2_UPSTREAM = "CACCGATT"
EMBED_BUILTIN_R2_SCHEME = PRIMER1 + "XX-" + (EMBED_R2_RESIDUAL + P7 + EMBED_R2_UPSTREAM)

# Internal substitutions (never at a read end) for the per-part max_errors tests.
P7_SUB = P7[:7] + ("C" if P7[7] != "C" else "A") + P7[8:]
SCAFF_R2_SUB = SCAFF_R2[:9] + ("A" if SCAFF_R2[9] != "A" else "C") + SCAFF_R2[10:]


# --- small helpers ----------------------------------------------------------

def _read_gz_records(path):
    """Return the FASTQ records in ``path`` (list of 4-line lists). Missing or
    empty files -> []."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with gzip.open(path, "rt") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    return [lines[i:i + 4] for i in range(0, len(lines), 4)]


def _cli_full(scheme, pairs, name_format, *extra, min_len=5):
    """Run the CLI over several read pairs; ``pairs`` is a list of
    ``(name, r1seq, r2seq)``.  Returns ``({R1:...,R2:...}, {R1:...,R2:...})``
    for the trimmed and discard outputs, handling empty/missing files."""
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/t1.fq", "w") as f1, open(f"{td}/t2.fq", "w") as f2:
            for name, r1seq, r2seq in pairs:
                f1.write(f"@{name}/1\n{r1seq}\n+\n{'I' * len(r1seq)}\n")
                f2.write(f"@{name}/2\n{r2seq}\n+\n{'I' * len(r2seq)}\n")
        p = subprocess.run(
            [CUTSEQ, "-A", scheme, "-O", f"{td}/o", "--rename", name_format,
             "-m", str(min_len), *extra, f"{td}/t1.fq", f"{td}/t2.fq"],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr

        def side(kind):
            return {s: _read_gz_records(f"{td}/o_{kind}_{s}.fastq.gz") for s in ("R1", "R2")}

        return side("trimmed"), side("discard")


def _run_yaml_pair(parts, r1, r2, name_format="{id}", q1=None, q2=None, settings=None):
    """Run the paired modifiers on a YAML part scheme (path handed to the
    public API); returns the two SequenceRecords."""
    import yaml
    with tempfile.TemporaryDirectory() as td:
        path = f"{td}/scheme.yaml"
        with open(path, "w") as fh:
            yaml.safe_dump({"parts": parts}, fh)
        return _run_paired(r1, r2, name_format=name_format, q1=q1, q2=q2,
                           scheme=path, settings=settings)


# --- the parent-flagged cases ------------------------------------------------

def test_full_readthrough_with_extra_bases_trimmed():
    """The p7 read-through adapter is followed by EXTRA non-adapter bases on both
    mates.  The trimming must consume the gate + extras + adjacent UMI/mask back
    to the insert (a suffix-only match would leave the trailing bases)."""
    r1 = XX2 + TRUN + INS + MASK8 + UMI8 + P7 + EXTRA_TAIL
    r2 = R2_UMI + _rc(MASK8) + RC_INS + _rc(TRUN) + _rc(XX2) + RC_P5 + EXTRA_TAIL
    r1o, r2o = _run_paired(r1, r2, name_format="{id}")
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence


def test_no_front_token_primer_insert_primer_full_readthrough_trimmed():
    """A primer:insert:primer scheme has NO read-start token (both outer
    adapters are recognised binding sites).  Full read-through must still trim
    the opposite arm — there must be no "progressed prerequisite" gate."""
    r1 = INS + P7
    r2 = RC_INS + _rc(PRIMER1)
    r1o, r2o = _run_paired(r1, r2, name_format="{id}", scheme=DSLIGATION_SCHEME)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence


def test_yaml_max_errors_far_site_substitution():
    """Per-part ``max_errors`` on the far (read-through) site: with 0.2 a single
    INTERNAL substitution is accepted and the UMI is captured; with 0 the same
    substitution is rejected, so the mirror cut does not fire while the read's
    own 5' cuts still happen."""
    parts = [
        {"type": "adapter", "seq": PRIMER1},
        {"type": "mask", "length": 2},
        {"type": "polytail", "base": "T"},
        {"type": "insert", "value": "-"},
        {"type": "mask", "length": 8},
        {"type": "capture", "length": 8, "label": "umi"},
        {"type": "adapter", "seq": P7, "max_errors": 0.2},
    ]
    r1 = XX2 + TRUN + INS + MASK8 + UMI8 + P7_SUB
    r2 = R2_UMI + _rc(MASK8) + RC_INS + _rc(TRUN) + _rc(XX2) + RC_P5

    # permissive: the one internal sub is tolerated -> full trim to insert; UMI captured
    r1o, r2o = _run_yaml_pair(parts, r1, r2, name_format="{id}_{r1.1}")
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence
    assert r1o.name == f"x/1_{UMI8}", r1o.name   # {r1.1} = top UMI captured

    # strict: same sub is NOT accepted -> mirror cut off; own 5' trims still fire
    strict = [
        ({**p, "max_errors": 0.0} if p.get("type") == "adapter" and p.get("seq") == P7 else dict(p))
        for p in parts
    ]
    r1s, r2s = _run_yaml_pair(strict, r1, r2, name_format="{id}_{r1.1}")
    assert r1s.sequence == INS + MASK8 + UMI8 + P7_SUB, r1s.sequence
    assert r1s.name == "x/1_", r1s.name            # {r1.1} empty: mirror UMI not captured
    assert r2s.sequence == RC_INS, r2s.sequence


def test_yaml_own_scaffold_sub_permissive_capture_correct():
    """An own-arm literal scaffold with one internal substitution is accepted
    under a permissive ``max_errors``; the adjacent UMI is captured correctly."""
    parts = [
        {"type": "adapter", "seq": PRIMER1},
        {"type": "mask", "length": 2},
        {"type": "insert", "value": "-"},
        {"type": "mask", "length": 4},
        {"type": "adapter", "seq": SCAFFOLD_TOP, "max_errors": 0.2},
        {"type": "capture", "length": 8, "label": "umi"},
        {"type": "adapter", "seq": P7},
    ]
    r1 = XX2 + INS
    r2 = R2_UMI + SCAFF_R2_SUB + MASK4_READ + RC_INS   # scaffold read with 1 internal sub
    r1o, r2o = _run_yaml_pair(parts, r1, r2, name_format="{id}_{1}")
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence       # permissive scaffold match -> X4 cut -> insert
    assert r2o.name == f"x/2_{R2_UMI}", r2o.name     # UMI captured correctly


def test_inline_after_umi_ensure_inline_barcode_rejects_mismatch():
    """A YAML inline barcode AFTER the UMI on the own arm: ``--ensure-inline-barcode``
    must recognise the matched exact object at the actual cursor (not the raw read
    start) and discard a pair whose inline region mismatches."""
    parts = [
        {"type": "adapter", "seq": PRIMER1},
        {"type": "mask", "length": 2},
        {"type": "capture", "length": 8, "label": "umi"},
        {"type": "inline", "seq": "ATCACG", "label": "bc"},
        {"type": "insert", "value": "-"},
        {"type": "mask", "length": 8},
        {"type": "capture", "length": 8, "label": "umi3"},
        {"type": "adapter", "seq": P7},
    ]
    good_r1 = XX2 + UMI8 + "ATCACG" + INS
    bad_r1 = XX2 + UMI8 + "TTTTTT" + INS                # inline region mismatches
    r2 = R2_UMI + _rc(MASK8) + RC_INS
    with tempfile.TemporaryDirectory() as td:
        path = f"{td}/scheme.yaml"
        with open(path, "w") as fh:
            import yaml
            yaml.safe_dump({"parts": parts}, fh)
        trimmed, discard = _cli_full(
            path, [("good", good_r1, r2), ("bad", bad_r1, r2)],
            "{id}_bc={bc}_umi={umi}_u3={umi3}", "--ensure-inline-barcode",
        )
    # the matching pair is kept, BOTH mates, naming the barcode (recognised at cursor)
    kept_name = f"bc=ATCACG_umi={UMI8}_u3={R2_UMI}"
    good_r1_recs, good_r2_recs = trimmed["R1"], trimmed["R2"]
    assert len(good_r1_recs) == 1, good_r1_recs
    assert len(good_r2_recs) == 1, good_r2_recs
    assert good_r1_recs[0][0] == f"@good/1_{kept_name}", good_r1_recs[0][0]
    assert good_r2_recs[0][0] == f"@good/2_{kept_name}", good_r2_recs[0][0]
    assert good_r1_recs[0][1] == INS, good_r1_recs[0][1]
    assert good_r2_recs[0][1] == RC_INS, good_r2_recs[0][1]
    # the mismatch pair is discarded, BOTH mates, out of trimmed and into discard
    assert all(not rec[0].startswith("@bad/") for rec in trimmed["R1"] + trimmed["R2"])
    assert all(not rec[0].startswith("@good/") for rec in discard["R1"] + discard["R2"])
    assert any(rec[0].startswith("@bad/1") for rec in discard["R1"]), discard["R1"]
    assert any(rec[0].startswith("@bad/2") for rec in discard["R2"]), discard["R2"]


def test_builtin_primer_repeated_distinct_boundaries_ambiguous():
    """The SAME builtin (known) primer at two distinct offsets in one arm is an
    ambiguous read-start boundary and must fail actionably (not only for custom
    primers)."""
    s = CutadaptConfig()
    s.auto_inline = False
    cs = _build_scheme(BUILTIN_AMBIGUOUS_SCHEME, s)
    with pytest.raises(ValueError) as excinfo:
        _scheme_modifiers(cs, paired=True, settings=s)
    assert "ambig" in str(excinfo.value).lower()


def test_embedded_builtin_r2_site_preserves_residual():
    """An embedded BUILTIN R2 binding site (P7) inside one merged uppercase token
    with a visible residual prefix (auto_inline off): the residual IS read on R2
    and consumed (a built-in site must not cause the whole token to be dropped),
    and the read-through arms are trimmed on both mates (no double-use of the
    upstream suffix beyond the site)."""
    s = CutadaptConfig()
    s.auto_inline = False
    s.r2_primer = _rc(P7)  # actual read-2 oligo for the builtin site

    # no read-through
    r1, r2 = _run_paired(XX2 + INS, _rc(EMBED_R2_RESIDUAL) + RC_INS,
                         name_format="{id}", scheme=EMBED_BUILTIN_R2_SCHEME, settings=s)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence  # residual consumed, not dropped

    # full read-through (upstream suffix included in the token)
    r1full = XX2 + INS + EMBED_R2_RESIDUAL + P7 + EMBED_R2_UPSTREAM
    r2full = _rc(EMBED_R2_RESIDUAL) + RC_INS + _rc(XX2) + _rc(PRIMER1)
    r1f, r2f = _run_paired(r1full, r2full, name_format="{id}",
                           scheme=EMBED_BUILTIN_R2_SCHEME, settings=s)
    assert r1f.sequence == INS, r1f.sequence
    assert r2f.sequence == RC_INS, r2f.sequence


def test_embedded_builtin_r2_site_auto_no_override():
    """Same embedded builtin R2 site, but with NO explicit ``--r2-primer``
    override: the builtin P7 site is recognised automatically (auto_inline still
    off), and the residual is consumed on both mates for no- and full
    read-through."""
    s = CutadaptConfig()
    s.auto_inline = False          # but no r2_primer override
    r1, r2 = _run_paired(XX2 + INS, _rc(EMBED_R2_RESIDUAL) + RC_INS,
                         name_format="{id}", scheme=EMBED_BUILTIN_R2_SCHEME, settings=s)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence

    r1full = XX2 + INS + EMBED_R2_RESIDUAL + P7 + EMBED_R2_UPSTREAM
    r2full = _rc(EMBED_R2_RESIDUAL) + RC_INS + _rc(XX2) + _rc(PRIMER1)
    r1f, r2f = _run_paired(r1full, r2full, name_format="{id}",
                           scheme=EMBED_BUILTIN_R2_SCHEME, settings=s)
    assert r1f.sequence == INS, r1f.sequence
    assert r2f.sequence == RC_INS, r2f.sequence


def test_single_end_left_primer_right_nonprimer_readthrough():
    """Single-end with a recognised left primer AND a right-hand NON-primer
    read-through scaffold: the read starts downstream of the left primer and the
    right read-through adapter is also trimmed."""
    s = CutadaptConfig()
    s.auto_inline = False
    out = _run_single(XX2 + INS + RIGHT_SCAFFOLD, name_format="{id}",
                      scheme=SE_RIGHT_SCHEME, settings=s)
    assert out.sequence == INS, out.sequence


def test_mixed_only_left_recognized_right_visible_scaffold():
    """Only the left side is a recognised primer; the right side is an EXPLICIT
    unknown (non-primer) read-visible scaffold.  R1 starts after the left binder
    and trims the right read-through; R2 begins at rc(right scaffold) (visible)
    and trims both that and the left read-through arm."""
    s = CutadaptConfig()
    s.auto_inline = False
    r1 = INS + RIGHT_SCAFFOLD
    r2 = _rc(RIGHT_SCAFFOLD) + RC_INS + _rc(PRIMER1)
    r1o, r2o = _run_paired(r1, r2, name_format="{id}", scheme=MIXED_SCHEME, settings=s)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence

# ============================================================================
# Multi-capture default naming (DECIDED: default concatenates ALL captured
# segments in SCHEME WRITTEN ORDER, not R2 execution order).  Same-R2 two N
# captures separated by a fixed known linker (type adapter), non-palindromic and
# distinct.  Source writer owns the default-renamer implementation; these encode
# the contract.
# ============================================================================

UMI_A_SEQ = "ACGTCGTA"
UMI_B_SEQ = "CATCGTAC"
LINKER_SEQ = "GTCAGTCA"
CELL_B_SEQ = "TGCATCGT"   # third (outer) capture, same R2 arm
CELL_L_SEQ = "ACGTTCGA"   # a capture on the LEFT (R1) arm
RC_A = _rc(UMI_A_SEQ)
RC_B = _rc(UMI_B_SEQ)
RC_LINKER = _rc(LINKER_SEQ)
RC_CELL_B = _rc(CELL_B_SEQ)


# --- 2-capture (single arm) + optional third / left capture -------------------

def _mc_parts(with_cell_b=False, with_cell_l=False):
    parts = [
        {"type": "adapter", "seq": PRIMER1},
    ]
    if with_cell_l:
        parts.append({"type": "capture", "length": 8, "label": "cell_l"})
    parts += [
        {"type": "mask", "length": 2},
        {"type": "insert", "value": "-"},
        {"type": "capture", "length": 8, "label": "umiA"},
        {"type": "adapter", "seq": LINKER_SEQ},
        {"type": "capture", "length": 8, "label": "umiB"},
    ]
    if with_cell_b:
        parts.append({"type": "capture", "length": 8, "label": "cell_bc"})
    parts.append({"type": "adapter", "seq": P7})
    return parts


def test_default_single_umi_unchanged():
    """The default for a single UMI is unchanged: ``{id}_{Umi}`` (here the
    contract scheme's single right UMI -> ``{id}_{rc(topUMI)}``)."""
    r1, r2 = _run_paired(_r1_nort(), _r2_nort(), name_format=None)
    assert r1.sequence == INS, r1.sequence
    assert r2.sequence == RC_INS, r2.sequence
    assert r1.name == f"x/1_{R2_UMI}", r1.name
    assert r2.name == f"x/2_{R2_UMI}", r2.name


def test_default_multicapture_same_r2_written_order():
    """Same-R2 two N captures (umIA + linker + umiB) read as rc(A) then linker
    then rc(A) on R2.  The DEFAULT name concatenates in SCHEME WRITTEN ORDER
    (pcA then umiB = rc(A)rc(B)), NOT the R2 execution order (rc(B)…rc(A))."""
    parts = _mc_parts()
    r1 = XX2 + INS
    r2 = RC_B + RC_LINKER + RC_A + RC_INS
    r1o, r2o = _run_yaml_pair(parts, r1, r2, name_format=None)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence
    assert r1o.name == f"x/1_{RC_A}{RC_B}", r1o.name
    assert r2o.name == f"x/2_{RC_A}{RC_B}", r2o.name


def test_explicit_grouped_multicapture_labels():
    """Explicit grouped template with a THIRD capture (same arm, distinct label):
    ``{id}_UMI:{1}{2}_BARCODE:{cell_bc}`` -> ``x/1_UMI:rc(A)rc(B)_BARCODE:rc(cellB)``
    on both mates; sequence + quality are the trimmed insert."""
    parts = _mc_parts(with_cell_b=True)
    r1 = XX2 + INS
    r2 = RC_CELL_B + RC_B + RC_LINKER + RC_A + RC_INS
    q1 = _qual(len(r1), seed=5)
    q2 = _qual(len(r2), seed=9)
    r1o, r2o = _run_yaml_pair(
        parts, r1, r2,
        name_format="{id}_UMI:{1}{2}_BARCODE:{cell_bc}",
        q1=q1, q2=q2,
    )
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence
    assert r1o.name == f"x/1_UMI:{RC_A}{RC_B}_BARCODE:{RC_CELL_B}", r1o.name
    assert r2o.name == f"x/2_UMI:{RC_A}{RC_B}_BARCODE:{RC_CELL_B}", r2o.name
    # qualities slice exactly to the insert only
    assert r1o.qualities == q1[len(XX2):], r1o.qualities
    assert r2o.qualities == q2[len(RC_CELL_B) + len(RC_B) + len(RC_LINKER) + len(RC_A):], r2o.qualities


def test_default_full_readthrough_multicapture_combines_once():
    """Full read-through: both mates observe the arm captures (R1 top-strand
    mirror and R2 rc).  The default must combine each capture ONCE (written
    order), not double-concatenate the mirrored values."""
    parts = _mc_parts()
    r1 = XX2 + INS + UMI_A_SEQ + LINKER_SEQ + UMI_B_SEQ + P7
    r2 = RC_B + RC_LINKER + RC_A + RC_INS + _rc(XX2) + _rc(PRIMER1)
    r1o, r2o = _run_yaml_pair(parts, r1, r2, name_format=None)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence
    assert r1o.name == f"x/1_{RC_A}{RC_B}", r1o.name   # not doubled
    assert r2o.name == f"x/2_{RC_A}{RC_B}", r2o.name


def test_default_multicapture_both_left_and_right_written_order():
    """Captures on BOTH arms: default concatenates in SCHEME WRITTEN ORDER
    (left cell_l first, then right umiA, umiB)."""
    parts = _mc_parts(with_cell_l=True)
    r1 = CELL_L_SEQ + XX2 + INS
    r2 = RC_B + RC_LINKER + RC_A + RC_INS
    r1o, r2o = _run_yaml_pair(parts, r1, r2, name_format=None)
    assert r1o.sequence == INS, r1o.sequence
    assert r2o.sequence == RC_INS, r2o.sequence
    assert r1o.name == f"x/1_{CELL_L_SEQ}{RC_A}{RC_B}", r1o.name
    assert r2o.name == f"x/2_{CELL_L_SEQ}{RC_A}{RC_B}", r2o.name


def test_duplicate_user_label_same_read_is_ambiguous():
    """Two captures on the SAME read/arm carrying the SAME user ``label`` is an
    ambiguous naming collision and must fail actionably (ValueError)."""
    parts = [
        {"type": "adapter", "seq": PRIMER1},
        {"type": "mask", "length": 2},
        {"type": "insert", "value": "-"},
        {"type": "capture", "length": 8, "label": "dup"},
        {"type": "adapter", "seq": LINKER_SEQ},
        {"type": "capture", "length": 8, "label": "dup"},
        {"type": "adapter", "seq": P7},
    ]
    with tempfile.TemporaryDirectory() as td:
        import yaml
        path = f"{td}/s.yaml"
        yaml.safe_dump({"parts": parts}, open(path, "w"))
        s = CutadaptConfig()
        with pytest.raises(ValueError):
            cs = _build_scheme(path, s)
            _scheme_modifiers(cs, paired=True, settings=s)


def test_autogenerated_barcode1_collision_with_explicit_barcode1():
    """An UNLABELED capture (auto ``barcode1``) next to an EXPLICIT ``barcode1``
    label is a collision and must fail actionably (reported, not silent)."""
    parts = [
        {"type": "adapter", "seq": PRIMER1},
        {"type": "mask", "length": 2},
        {"type": "insert", "value": "-"},
        {"type": "capture", "length": 8},                    # auto -> barcode1
        {"type": "adapter", "seq": LINKER_SEQ},
        {"type": "capture", "length": 8, "label": "barcode1"},  # collision
        {"type": "adapter", "seq": P7},
    ]
    with tempfile.TemporaryDirectory() as td:
        import yaml
        path = f"{td}/s.yaml"
        yaml.safe_dump({"parts": parts}, open(path, "w"))
        s = CutadaptConfig()
        with pytest.raises(ValueError):
            cs = _build_scheme(path, s)
            _scheme_modifiers(cs, paired=True, settings=s)


def test_incomplete_inner_fragment_other_complete_explicit_refs():
    """R2 is truncated inside the INNER (umiA) capture: the explicit reference
    for the incomplete inner capture is EMPTY while the complete outer (umiB)
    keeps its value, using separators: ``{id}_inner:{1}_outer:{2}``.  (No
    claim about the not-yet-finalised default incomplete policy.)"""
    parts = _mc_parts()
    r1 = XX2 + INS
    r2 = RC_B + RC_LINKER + RC_A[:3]      # inner rc(A) present only 3 bp
    r1o, r2o = _run_yaml_pair(parts, r1, r2, name_format="{id}_inner:{1}_outer:{2}")
    assert r1o.sequence == INS, r1o.sequence
    assert r1o.name == f"x/1_inner:_outer:{RC_B}", r1o.name
    assert r2o.name == f"x/2_inner:_outer:{RC_B}", r2o.name
