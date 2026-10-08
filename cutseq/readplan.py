#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Read-plan routing for primer-boundary library schemes.

cutseq's legacy ``compile_tokens`` treats every scheme token as read-visible on
BOTH reads (written-side at the 5' end, mirrored at the other read's 3' end).
That is exactly right for ``SCG``/``DBiT`` style schemes where the outer
adapters are in-read scaffold, but it is WRONG for a library whose outer
adapters are genuine *sequencing-primer binding sites*: the sequencing primer
anneals UPSTREAM of the read, so the read does not contain the binder at all — it
starts just inside it.

This module builds that read plan.  The model is:

* **Boundary** — a recognised sequencing-site primer locates, per READ side, a
  strand-specific binding site and the read's 3' extension boundary *in the
  scheme* (not in the FASTQ).  A ``Boundary`` records ``(side, token_index,
  within-token offset, strand, role, site_seq, primer_seq)``.  Everything at or
  after the boundary offset on the read's own arm is read-visible (residual
  scaffold inside the same token is preserved); everything strictly before it is
  upstream and skipped on THAT read — but a token upstream of the boundary on the
  OPPOSITE read's read-through walk is never processed as an inward cut.
* **Per-read plan** — for each of R1 and R2 independently, an ordered list of
  steps built from the recognised boundary (or a read-visible fallback when no
  primer was found).  A read plan is NOT gated on both primers being present:
  one recognised side coexists with a read-visible other side, and single-end
  with a recognised R1 site also starts downstream.
* **Arm-local executor** — a picklable, per-read ``_PlanExecutor`` owns the
  ordered steps for one read and a *local* cursor.  Each read-through step is
  gated on the immediately-preceding step having produced a match at the current
  cursor (progression), never on a cumulative ``info.matches`` scan and never on
  ``len >= 50``.  Capture observations are recorded with provenance
  (start/end/as-read/status) so an incomplete capture is never emitted as if it
  were a complete one.

Primer identification
---------------------
The built-in ``SEQUENCING_SITES`` map (see ``cutseq/primers``) is direction- and
role-aware, plus explicitly supplied ``--r1-primer`` / ``--r2-primer`` (the
actual read oligos) override the DB.  A recognised site must map to a genuine
3' extension boundary: for the R1 (left) site only the primer's 3' end defines
where the read starts; for the R2 (right) site only the top-strand site's prefix
defines it.  An arbitrary primer prefix is NOT a boundary.  When several distinct
plausible boundaries resolve for the same side the compiler fails actionably
rather than guessing; same-boundary aliases merge.

Short-capture completeness policy
---------------------------------
An incomplete read (shorter than a declared capture/mask) records exactly the
bases present and never fabricates the missing tail.  A read-START capture
records ``min(len(read), n)`` bases and consumes the whole read (the empty read
is then discarded by ``--min-length``); the other mate is untouched.  A
read-THROUGH capture only fires once the read actually read through, and
records/trims the contiguous present bases.  No automatic canonicalisation is
applied (``rc()``/``rev()``/``comp()`` stay explicit).

If a meaningful, unresolved ambiguity prevents locating a single boundary the
compiler raises an actionable error instead of guessing.  Schemes with no
recognised primer fall back to the legacy read-visible per-token emission.
"""

import cutadapt.modifiers as _mods


def _g():
    import cutseq.grammar as gr
    return gr


def _rc(seq):
    return _g()._rc(seq)


def _norm(v):
    return str(v).upper().replace("U", "T").translate(
        str.maketrans("", "", " *:+-.'\u2032\u00b0")
    )


# --- boundary -------------------------------------------------------------


class Boundary:
    """Where one read starts in the scheme, given a recognised sequencing site.

    ``token_index`` / ``token_offset`` locate the read start in the written arm:
    the read's own-arm content begins at ``token_offset`` within
    ``arm[token_index]`` (and continues through the rest of that token and the
    following tokens).  ``strand`` is ``'top'`` for R1 (reads the top strand
    5'->3') or ``'bot'`` for R2 (reads the bottom strand from the outer end
    inward).  ``role`` is ``'R1'`` or ``'R2'``; ``site_seq`` is the recognised
    top-strand site; ``primer_seq`` is the read oligo (for R2 this is rc of the
    top-strand site).
    """

    __slots__ = (
        "side", "token_index", "token_offset", "strand", "role",
        "site_seq", "primer_seq",
    )

    def __init__(self, side, token_index, token_offset, strand, role,
                 site_seq, primer_seq):
        self.side = side
        self.token_index = token_index
        self.token_offset = token_offset
        self.strand = strand
        self.role = role
        self.site_seq = site_seq
        self.primer_seq = primer_seq

    def __repr__(self):
        return (f"Boundary({self.side}, tok{self.token_index}+"
                f"{self.token_offset}, {self.strand}, {self.role})")


def _resolve_boundary(arm, side, custom=None):
    """Resolve the recognised sequencing-site boundary on ``arm``.

    Returns ``(Boundary, upstream_tokens)`` or ``None`` when no site maps to a
    genuine extension boundary.  When a CUSTOM primer is supplied it OVERRIDES
    the built-in DB: it must locate a single distinct read-start offset, else
    the resolve fails actionably (absent -> ValueError; ambiguous -> ValueError
    containing ``ambig``).  Same-boundary aliases (distinct site records giving
    the same boundary offset) merge; distinct boundaries raise.
    """
    from .primers import (_rc as _prc, _norm as _pn, SEQUENCING_SITES,
                          _site_boundaries)

    # --- explicit custom primers override the DB ---------------------------
    if custom is not None:
        return _resolve_custom(arm, side, custom)

    idxs = [i for i, t in enumerate(arm) if t.kind == "adp"]
    # Collect EVERY boundary candidate across the whole arm (including sites
    # separated by N/X tokens and sites embedded within a single merged token),
    # not just the first token that happens to match.
    candid = []  # (token_index, offset, site)
    for i in idxs:
        value = _pn(arm[i].value)
        for _name, (role, site) in SEQUENCING_SITES.items():
            if role != side:
                continue
            for off in _site_boundaries(value, site, side):
                candid.append((i, off, site))
    if not candid:
        return None
    # Same-coordinate aliases (distinct site records -> the SAME (token,offset))
    # merge; keep the longest site at each coordinate.  >1 distinct coordinate is
    # genuinely ambiguous (an incompatible multi-assay map) and fails actionably.
    best = {}
    for i, off, site in candid:
        k = (i, off)
        if k not in best or len(site) > len(best[k]):
            best[k] = site
    if len(best) > 1:
        raise ValueError(
            "cutseq: multiple distinct " + ("R1" if side == "R1" else "R2") +
            " primer boundaries (ambiguous) on the " +
            ("left" if side == "R1" else "right") +
            " arm; specify --r1-primer / --r2-primer to disambiguate"
        )
    (ti, off), site = next(iter(best.items()))
    value = _pn(arm[ti].value)
    # The read-through GATE is the ACTUAL observed site fragment inside the token
    # (e.g. the short read-adjacent adapter actually written), not the longest DB
    # alias that happens to share the same boundary offset.
    if side == "R1":
        observed = value[max(0, off - len(site)):off]
        return (Boundary("R1", ti, off, "top", "R1", observed, observed),
                arm[:ti])
    observed = value[off:off + len(site)]
    return (Boundary("R2", ti, off, "bot", "R2", observed, _prc(observed)),
            arm[ti + 1 :])


def _resolve_custom(arm, side, custom):
    """Resolve a boundary for an explicitly supplied custom read-oligo.

    ``custom`` is the actual read primer 5'->3'.  For ``side == 'R1'`` the read
    starts at the primer's 3' extension end; for ``'R2'`` the top-strand site is
    ``rc(custom)`` and the read starts at its left (insert-adjacent) edge.  All
    occurrences are collected across the whole arm; 0 -> absent, >1 distinct
    offsets -> ambiguous (both raise actionably).
    """
    from .primers import _rc as _prc, _norm as _pn, MIN_PRIMER_MATCH

    primer = _pn(custom)
    # A custom primer must locate a genuine boundary.  A short custom primer is
    # accepted verbatim ("full shorter user custom"); a longer one locates a
    # boundary only via a full occurrence or a terminal stem of at least the
    # curated fragment strength (don't credit a 1 bp fragment of a long primer).
    min_frag = min(len(primer), MIN_PRIMER_MATCH)
    hits = []  # (token_index, token_offset, site_seq)
    for i, t in enumerate(arm):
        if t.kind != "adp":
            continue
        value = _pn(t.value)
        if side == "R1":
            # boundary at the primer 3' end wherever the primer appears fully
            # (a +1 scan catches overlapping occurrences).
            pos = value.find(primer)
            while pos != -1:
                hits.append((i, pos + len(primer), primer))
                pos = value.find(primer, pos + 1)
            # a terminal stem shorter than the primer (>= min_frag)
            if len(value) < len(primer) and primer.endswith(value) and \
                    len(value) >= min_frag:
                hits.append((i, len(value), primer))
        else:
            site = _prc(primer)
            pos = value.find(site)
            while pos != -1:
                hits.append((i, pos, site))
                pos = value.find(site, pos + 1)
            if len(value) < len(site) and site.startswith(value) and \
                    len(value) >= min_frag:
                hits.append((i, 0, site))
    if not hits:
        raise ValueError(
            f"cutseq: supplied {'--r1-primer' if side == 'R1' else '--r2-primer'} "
            f"{custom!r} was not found in the scheme"
        )
    # deduplicate distinct boundaries
    distinct = {}
    for ti, off, site in hits:
        distinct[(ti, off)] = site
    if len(distinct) > 1:
        raise ValueError(
            "cutseq: custom primer resolves to multiple distinct read-start "
            f"boundaries on the {'left' if side == 'R1' else 'right'} arm "
            "(ambiguous); it should occur once"
        )
    (ti, off), site = next(iter(distinct.items()))
    if side == "R1":
        return (Boundary("R1", ti, off, "top", "R1", site, primer), arm[:ti])
    return (Boundary("R2", ti, off, "bot", "R2", site, primer), arm[ti + 1 :])



# --- read-step model ------------------------------------------------------


class _Step:
    """One ordered action on one read.

    ``kind``: ``'start'`` (read-start positional trim/capture/adapter) or
    ``'through'`` (read-through, gated on progression).  ``act`` is ``'trim'``,
    ``'capture'``, ``'poly5'``, ``'poly3'``, ``'adapter'`` or ``'inline'``.
    ``length`` is a fixed base count (mask/capture), 0 for variable runs/adapters.
    ``label`` is the capture id.  ``seq`` is a fixed adapter/inline sequence.
    ``opts`` is the per-token match options (never silently dropped).
    """

    __slots__ = ("kind", "act", "length", "label", "seq", "opts",
                 "capture_id")

    def __init__(self, kind, act, length=0, label=None, seq=None, opts=None):
        self.kind = kind
        self.act = act
        self.length = length
        self.label = label
        self.seq = seq
        self.opts = opts or {}
        self.capture_id = label

    def __repr__(self):
        base = f"ReadStep({self.kind}/{self.act}"
        if self.length:
            base += f",len={self.length}"
        if self.label:
            base += f",name={self.label}"
        return base + ")"


class ReadPlan:
    """Explicit per-read routing plan for a primer-boundary scheme.

    ``r1_boundary``/``r2_boundary`` are ``Boundary`` or ``None`` (None == that
    side is read-visible).  ``r1_steps``/``r2_steps`` are the ordered per-read
    steps.  ``r1_upstream``/``r2_upstream`` are the tokens skipped on that read
    (upstream of its boundary).  ``mode`` is ``'readplan'`` when at least one
    side had a recognised primer, else ``'read-visible'``.  ``r1_is_plan``/
    ``r2_is_plan`` record which side was routed through the read-plan engine so
    the caller can mix in the legacy per-token emission for the other side.
    """

    __slots__ = (
        "r1_boundary", "r2_boundary", "r1_steps", "r2_steps",
        "r1_capture_ids", "r2_capture_ids", "mode", "r1_upstream",
        "r2_upstream", "r1_is_plan", "r2_is_plan",
        "r1_own_arm_adapters", "r2_own_arm_adapters",
    )

    def __init__(self, r1_boundary, r2_boundary, r1_steps, r2_steps,
                 r1_upstream, r2_upstream):
        self.r1_boundary = r1_boundary
        self.r2_boundary = r2_boundary
        self.r1_steps = r1_steps
        self.r2_steps = r2_steps
        self.r1_upstream = r1_upstream
        self.r2_upstream = r2_upstream
        # The exact own-arm (read-start) adapter objects each read's executor
        # matches with, so ``r2_arm_adapters`` / ``NoCassette`` can use identity.
        self.r1_own_arm_adapters = []
        self.r2_own_arm_adapters = []
        self.r1_capture_ids = [s.label for s in r1_steps if s.label]
        self.r2_capture_ids = [s.label for s in r2_steps if s.label]
        self.mode = (
            "readplan"
            if (r1_boundary is not None or r2_boundary is not None)
            else "read-visible"
        )
        # Whenever a plan exists BOTH reads are routed through an arm-local
        # executor (the unrecognised side uses its read-visible own tokens), so
        # both are flagged "planned" for the caller.
        self.r1_is_plan = self.mode == "readplan"
        self.r2_is_plan = self.mode == "readplan"

    def __repr__(self):
        return (f"ReadPlan(mode={self.mode}, r1b={self.r1_boundary},"
                f" r2b={self.r2_boundary})")


# --- arm-local executor ---------------------------------------------------


class _Observation:
    """Provenance of one recorded capture on a read."""

    __slots__ = ("label", "start", "end", "seq", "complete", "side")

    def __init__(self, label, start, end, seq, complete, side):
        self.label = label
        self.start = start
        self.end = end
        self.seq = seq
        self.complete = complete
        self.side = side

    def __repr__(self):
        c = "ok" if self.complete else "short"
        return f"Observation({self.label}[{self.start}:{self.end}]={self.seq!r},{c})"


class _PlanExecutor(_mods.SingleEndModifier):
    """Per-read executor that walks one read's plan with a local cursor.

    A single modifier per read; keeps ONLY local state for the current read
    (a local cursor, observations and recorded matches) and never writes mutable
    cursor state onto ``self``, so an executor is safe to share across threads /
    workers.  Read-through steps are gated on the immediately-preceding read-\nthrough result (a far adapter matched) -- never on ``info.matches`` and never
    on ``len >= 50``.  Matching uses real cutadapt adapters (``PrefixAdapter``
    for own-arm scaffolds at the 5' cursor, ``BackAdapter`` for the first far
    read-through site, ``SuffixAdapter`` for subsequent inward fixed linkers) so
    per-token ``max_errors`` / ``indels`` / wildcards and overlap are honoured,
    and the actual matched adapter (identity) is recorded for the filters
    (``--ensure-inline-barcode``, ``--require-cassette``).
    """

    def __init__(self, steps, side, force_anywhere=False):
        self.steps = list(steps)
        self.side = side
        self.force_anywhere = force_anywhere
        # Each step gets its OWN adapter object (honouring its per-token options)
        # so two same-sequence inline steps with different options are NOT
        # conflated, and the exposed ``adapters`` are exactly what we match with.
        self._adps, self._inline_adapters = self._build_adapters()
        self._is_inline = bool(self._inline_adapters)

    @property
    def adapters(self):
        """The read-visible inline-barcode adapter objects this read plan matches
        (used by ``_collect_inline`` / ``inline_adapters`` / IsUntrimmedAny)."""
        return list(self._inline_adapters)

    def _build_adapters(self):
        """Return ``(adapter_by_step_index, inline_adapters)``.

        * read-start ``adapter``/``inline`` -> ``PrefixAdapter`` (anchored 5').
        * the FIRST read-through ``adapter``/``inline`` -> ``BackAdapter`` (the
          opposite-arm sequencing site, tolerant of extra bases after it).
        * any LATER read-through ``adapter``/``inline`` -> ``SuffixAdapter`` (a
          fixed linker already at the advanced 3' cursor).
        Inline read-start adapters are also exposed for ``--ensure-inline-barcode``,
        each built from its own token options so identity is step-specific.
        """
        from cutadapt.adapters import BackAdapter, PrefixAdapter, SuffixAdapter

        out = {}
        inline = []
        through_adapter_seen = False
        for idx, step in enumerate(self.steps):
            if step.act not in ("adapter", "inline"):
                out[idx] = None
                continue
            if step.kind == "start":
                adp = PrefixAdapter(
                    step.seq,
                    **self._adapter_kw(step, min_overlap=4),
                )
                out[idx] = adp
                # Only read-START (own-arm) inline barcodes are exposed for
                # ``--ensure-inline-barcode``; a read-through inline is a gate,
                # not an expected 5' barcode, so it must not require a match.
                if step.act == "inline":
                    inline.append(adp)
            elif step.kind == "through":
                if not through_adapter_seen:
                    adp = BackAdapter(
                        step.seq,
                        **self._adapter_kw(step, min_overlap=10, back=True),
                    )
                    through_adapter_seen = True
                else:
                    adp = SuffixAdapter(
                        step.seq,
                        **self._adapter_kw(step, min_overlap=4),
                    )
                out[idx] = adp
        return out, inline

    def own_arm_adapters(self):
        """The own-arm (read-start) adapter/inline objects this read matches."""
        return [
            a for idx, a in sorted(self._adps.items())
            if a is not None and self.steps[idx].kind == "start"
        ]

    def step_adapters(self):
        """``[(step, adapter)]`` for every adapter/inline step (own-arm + read-through)."""
        out = []
        for idx, step in enumerate(self.steps):
            a = self._adps.get(idx)
            if a is not None:
                out.append((step, a))
        return out

    def _adapter_kw(self, step, min_overlap=4, back=False):
        """Cutadapt adapter kwargs from a step's per-token options + defaults.

        An explicit ``False`` for wildcards / indels is preserved (not replaced
        by the cutadapt default); ``force_anywhere`` is honoured when set on the
        token or globally.
        """
        o = step.opts or {}
        kw = {
            "max_errors": o.get("max_errors", 0.2),
            "min_overlap": o.get("min_overlap", min_overlap),
        }
        if "read_wildcards" in o:
            kw["read_wildcards"] = o["read_wildcards"]
        if "adapter_wildcards" in o:
            kw["adapter_wildcards"] = o["adapter_wildcards"]
        if "indels" in o:
            kw["indels"] = o["indels"]
        if o.get("force_anywhere") or self.force_anywhere:
            kw["force_anywhere"] = True
        return kw

    def _run(self, read, info):
        s = read.sequence
        n = len(s)
        start = 0          # consumed from the 5' (read-start) end
        end = n            # consumed from the 3' (read-through) end
        obs = []
        through_anchor = False      # a read-through adapter matched this walk
        block_through = False       # the far gate failed: no re-anchoring
        scaffold_failed = False     # an own-arm scaffold failed to match
        for idx, step in enumerate(self.steps):
            if step.kind == "through":
                if block_through:
                    continue
                if step.act in ("adapter", "inline"):
                    if not through_anchor:
                        # first read-through site: match anywhere in the current
                        # [start:end] window (not the consumed 5' prefix)
                        m = self._match_back(step, s, start, end, idx)
                        if m is None:
                            block_through = True
                            continue
                        # A partial match that misses INWARD adapter bases
                        # (``astart > 0``) means the read did not reach the site's
                        # inner edge, so we must NOT infer the UMI/mask boundary:
                        # trim the matched outer part but disable dependent cuts.
                        end = max(start, min(end, start + m.rstart))
                        through_anchor = m.astart == 0
                        self._record_match(info, m)
                        if step.label:
                            obs.append(_Observation(
                                step.label, start + m.rstart, start + m.rstop,
                                s[start + m.rstart:start + m.rstop], True, self.side))
                        continue
                    # a later inward fixed linker: anchored at the new 3' cursor
                    m = self._match_suffix(step, s, start, end, idx)
                    if m is None:
                        block_through = True
                        continue
                    end = max(start, min(end, start + m.rstart))
                    through_anchor = through_anchor and m.astart == 0
                    self._record_match(info, m)
                    if step.label:
                        obs.append(_Observation(
                            step.label, start + m.rstart, start + m.rstop,
                            s[start + m.rstart:start + m.rstop], True, self.side))
                    continue
                if not through_anchor:
                    continue
                if step.act == "capture":
                    n_ = step.length
                    e = end
                    st = e - n_ if e - n_ >= start else start
                    cap = s[st:e]
                    obs.append(_Observation(step.label, st, e, cap,
                                            e - st == n_, self.side))
                    end = st
                elif step.act == "trim":
                    n_ = step.length
                    e = end
                    st = e - n_ if e - n_ >= start else start
                    end = st
                elif step.act == "poly3":
                    win = s[start:end]
                    cut = _g()._poly_tail_trim_index(win, step.seq)
                    # A too-short 3' run (< min_len / min_overlap) must NOT be
                    # consumed: honour the explicit floor, else the native
                    # score-cutoff applies.  ``cut`` is the index where the 3'
                    # run begins, so the run length is ``len(win) - cut``.
                    floor = _poly_run_min(step.opts)
                    if len(win) - cut >= floor:
                        end = max(start, min(end, start + cut))
            else:  # read-start
                # Once an own-arm scaffold fails to match, EVERY subsequent
                # own-arm step stays blocked (no re-anchoring / re-enabling).
                if scaffold_failed:
                    continue
                if step.act in ("adapter", "inline"):
                    m = self._match_prefix_step(step, s, start, idx)
                    if m is None:
                        # failed read-visible scaffold: block downstream
                        # positional cuts (they must not consume insert bases)
                        scaffold_failed = True
                        continue
                    end = n  # unchanged
                    # Record the actual match so the filter sees it even when the
                    # matched region is not at the raw read start.
                    self._record_match(info, m)
                    if step.label:
                        obs.append(_Observation(
                            step.label, start + m.rstart, start + m.rstop,
                            s[start + m.rstart:start + m.rstop], True, self.side))
                    start += m.rstop
                elif step.act == "trim":
                    start += min(step.length, max(0, n - start))
                elif step.act == "capture":
                    n_ = step.length
                    avail = n - start
                    if avail < n_:
                        cap, complete, e = s[start:], False, n
                    else:
                        cap, complete, e = s[start : start + n_], True, start + n_
                    obs.append(_Observation(step.label, start, e, cap,
                                            complete, self.side))
                    start = e
                elif step.act == "poly5":
                    cut = _g()._poly_head_trim_index(s[start:], step.seq)
                    # A too-short 5' run (< min_len / min_overlap) must NOT be
                    # consumed: honour the explicit floor.  ``cut`` IS the number
                    # of leading run bases, so compare it directly.
                    if cut >= _poly_run_min(step.opts):
                        start += cut
        # ``s[start:end]`` is naturally empty when start > end; never resurrect
        # the whole (untouched) read for an over-consumed position.
        final = s[start:end] if start <= end else ""
        return final, obs, start, end

    def _match_prefix_step(self, step, s, start, idx):
        """Match an own-arm adapter/inline at the 5' cursor; Match or None."""
        adp = self._adps[idx]
        if adp is None:
            return None
        return adp.match_to(s[start:])

    def _match_back(self, step, s, start, end, idx):
        """Match the first read-through far site in the current ``[start:end]``
        window (reads through the whole consumed 5' prefix are NOT re-scanned)."""
        adp = self._adps[idx]
        if adp is None:
            return None
        return adp.match_to(s[start:end])

    def _match_suffix(self, step, s, start, end, idx):
        """Match an inward fixed linker anchored at the current 3' cursor."""
        adp = self._adps[idx]
        if adp is None:
            return None
        return adp.match_to(s[start:end])

    def _record_match(self, info, m):
        """Append a real matched adapter so the filters/stats see it (identity)."""
        if m is not None:
            info.matches.append(m)

    def __call__(self, read, info):
        final, obs, start, end = self._run(read, info)
        for o in obs:
            # Only COMPLETE captures are advertised under an ordinary ``{n}`` /
            # label reference.  An incomplete capture (read shorter than the
            # declared capture) is kept as a diagnostic observation but is NOT
            # recorded as if it were a valid full-length one, so the rename
            # resolves it to empty.  The read itself is still consumed.
            if o.complete:
                _g()._record_capture(info, o.label, o.seq)
                if o.start == 0:
                    info.cut_prefix = o.seq
        if start > end:
            start, end = 0, 0
        if final == read.sequence:
            return read
        from dnaio import SequenceRecord

        # Cut the quality string with the SAME start/end offsets as the sequence
        # so quality slicing stays in sync with the read-position trims.
        q = read.qualities
        sl = q[start:end] if q else None
        return SequenceRecord(read.name, final, sl)

    def __repr__(self):
        # Expose the ordered steps (kinds / actions / fixed lengths / capture
        # labels) so a dry-run or summary shows the read plan actually applied.
        side = "R1" if self.side == 0 else "R2"
        parts = []
        for s in self.steps:
            tag = s.kind.replace("start", "5'").replace("through", "3'")
            d = f"{tag}:{s.act}"
            if s.length:
                d += f"({s.length})"
            if s.label:
                d += f",name={s.label}"
            parts.append(d)
        return f"PlanExecutor({side}, steps={len(self.steps)} [{' '.join(parts)}])"


# --- step construction ----------------------------------------------------


def _token_opts(tok):
    """Per-token match options honoured by the executor (never dropped).

    Includes ``force_anywhere`` and preserves an explicit ``False`` for the
    wildcard / indel flags so a ``False`` written in a YAML part is enforced
    (not silently replaced by the cutadapt default).
    """
    return {
        k: tok.options[k]
        for k in ("max_errors", "min_overlap", "min_len", "read_wildcards",
                  "adapter_wildcards", "indels", "force_anywhere")
        if k in tok.options
    }


def _poly_run_min(opts):
    """The minimum homopolymer run length a poly-5'/poly-3' step may trim.

    Per-token ``min_len`` and ``min_overlap`` both act as a floor on the run a
    poly-5'/poly-3' step will consume; with neither set the native
    ``poly_a_trim_index`` threshold governs (a run must be non-trivial to be
    trimmed at all), so an unconfigured poly is left unchanged.  Returns an
    int floor (0 = no explicit floor -> honour the native score-cutoff only).
    """
    o = opts or {}
    floor = 0
    for k in ("min_len", "min_overlap"):
        val = o.get(k)
        if val is None:
            continue
        try:
            floor = max(floor, int(float(val)))
        except (TypeError, ValueError):
            continue
    return floor


def _poly_check_opts(tok):
    """Validate a poly token's options before compiling it into a read step.

    The read-plan poly-5'/poly-3' trims reuse cutadapt's *score-based*
    ``poly_a_trim_index``, which is NOT an exact-match adapter and therefore
    cannot honour an arbitrary per-token ``max_errors``.  An explicitly
    supplied non-default ``max_errors`` on a poly token is silently ignored by
    that native routine, so rather than claim it is honoured we fail actionably.
    The default ``0.2`` (or an absent key) is a no-op and is left alone.
    """
    me = (tok.options or {}).get("max_errors")
    if me is not None:
        try:
            me_f = float(me)
        except (TypeError, ValueError):
            me_f = None
        if me_f is not None and me_f != 0.2:
            raise ValueError(
                "cutseq: poly tail/repeat options do not support arbitrary "
                f"'max_errors' ({me!r}); the native homopolymer trimmer is "
                "score-based. Set 'min_len' / 'min_overlap' to control how long "
                "a run is trimmed, or remove 'max_errors'."
            )


def _start_step(t, r2=False):
    """Read-start step for a token.  ``r2`` rc's the adapter/inline/poly bases so
    the executor matches them against the bottom-strand read."""
    def v(seq):
        return _rc(seq) if r2 else seq
    if t.kind == "mask":
        return _Step("start", "trim", length=t.value)
    if t.kind == "capture":
        return _Step("start", "capture", length=t.value, label=t.label)
    if t.kind == "polytail":
        _poly_check_opts(t)
        return _Step("start", "poly5", seq=v(t.value), opts=_token_opts(t))
    if t.kind in ("adp", "back"):
        return _Step("start", "adapter", seq=v(t.value), opts=_token_opts(t))
    if t.kind == "inline":
        return _Step("start", "inline", seq=v(t.value), label=t.label,
                     opts=_token_opts(t))
    return _Step("start", "trim", length=0)


def _through_step(t, r2=False):
    """Read-through step for a token (``r2`` rc's sequence bases)."""
    def v(seq):
        return _rc(seq) if r2 else seq
    if t.kind == "mask":
        return _Step("through", "trim", length=t.value)
    if t.kind == "capture":
        return _Step("through", "capture", length=t.value, label=t.label)
    if t.kind == "polytail":
        _poly_check_opts(t)
        return _Step("through", "poly3", seq=v(t.value), opts=_token_opts(t))
    if t.kind in ("adp", "back"):
        return _Step("through", "adapter", seq=v(t.value), opts=_token_opts(t))
    if t.kind == "inline":
        return _Step("through", "inline", seq=v(t.value), label=t.label,
                     opts=_token_opts(t))
    return _Step("through", "trim", length=0)


def _read_start_steps(boundary, arm):
    """Read-start steps for a boundary: within-token residual + following tokens.

    R1 reads the top strand in written order; R2 reads the bottom strand from the
    outer end inward (reversed written order, rc'd).  The residual inside the
    boundary token (the part of an adapter token at/after ``token_offset`` on
    R1, or before the site on R2) is preserved and trimmed as an anchored mate.
    """
    side = boundary.side
    r2 = side == "R2"
    i, off = boundary.token_index, boundary.token_offset
    steps = []
    if not r2:
        btok = arm[i]
        if btok.kind == "adp" and off < len(btok.value):
            steps.append(_Step("start", "adapter", seq=btok.value[off:],
                               opts=_token_opts(btok)))
        for t in arm[i + 1 :]:
            steps.append(_start_step(t))
    else:
        if off > 0:
            btok = arm[i]
            steps.append(_Step("start", "adapter", seq=_rc(btok.value[:off]),
                               opts=_token_opts(btok)))
        for t in reversed(arm[:i]):
            steps.append(_start_step(t, r2=True))
    return steps


def _read_through_steps(boundary, arm, read_side):
    """Read-through steps for ``read_side`` walking the OPPOSITE arm, from the far
    boundary site inward (tokens upstream of the far site are excluded).

    ``boundary`` is the FAR side's boundary, or ``None``.  The boundary's own
    token is SPLIT so the read-through GATE is the observed primer site region
    only (not the whole merged token): the inward residual is a separate step and
    the upstream remainder is trimmed with the gate (never required for a match),
    so a partially-sequenced site neither mis-gates nor cuts the residual twice.
    """
    if boundary is None:
        if read_side == "R1":
            return [_through_step(t) for t in reversed(arm)]
        return [_through_step(t, r2=True) for t in arm]
    i = boundary.token_index
    btok = arm[i]
    site = boundary.site_seq
    off = boundary.token_offset
    val = btok.value
    # The read-through gate is the SITE region (rc'd for R2).
    if read_side == "R1":
        # boundary token (right arm, top strand): residual + site + upstream.
        steps = []
        gate = _Step("through", "adapter", seq=site, opts=_token_opts(btok))
        steps.append(gate)
        if off > 0:
            steps.append(_Step("through", "adapter", seq=val[:off],
                               opts=_token_opts(btok)))
        # tokens between the boundary and the insert (index i..0 in written), then
        # the insert-adjacent remainder of the arm (outermost-first).
        for t in reversed(arm[:i]):
            steps.append(_through_step(t))
        return steps
    # R2: left arm bottom-strand, written order rc'd.  Boundary token = upstream +
    # site + residual; the gate is rc(site); the inward residual is rc(val[off:]).
    steps = []
    gate = _Step("through", "adapter", seq=_rc(site), opts=_token_opts(btok))
    steps.append(gate)
    if off < len(val):
        steps.append(_Step("through", "adapter", seq=_rc(val[off:]),
                           opts=_token_opts(btok)))
    for t in arm[i + 1 :]:
        steps.append(_through_step(t, r2=True))
    return steps


def _build_plan(left, right, r1b, r2b, orientation=None):
    """Construct a coherent per-read `ReadPlan` once ANY side resolves.

    ``r1b``/``r2b`` are ``(Boundary, upstream)`` or ``None``.  Returns ``None``
    when NEITHER side resolved (fully read-visible -> legacy per-token path).
    When at least one side is a recognised primer, BOTH reads are routed through
    an arm-local plan: the recognised side starts at its binding boundary, the
    unrecognised side keeps ALL of its own read-visible tokens (from the outer
    end inward) and gets the opposite correct binder's read-through -- so one
    recognised side never silently collapses the other mate to a legacy chain
    that gates/cuts incorrectly.
    """
    if r1b is None and r2b is None:
        return None
    r1bnd = r1b[0] if r1b else None
    r2bnd = r2b[0] if r2b else None
    r1_up = r1b[1] if r1b else []
    r2_up = r2b[1] if r2b else []
    # R1: start at its boundary (or read-visible whole left arm) + right read-through.
    if r1bnd is not None:
        r1_steps = _read_start_steps(r1bnd, left)
    else:
        r1_steps = _read_visible_start_steps(left, r2=False)
    r1_steps += _read_through_steps(r2bnd, right, "R1")
    # R2: the mirror of the above.
    if r2bnd is not None:
        r2_steps = _read_start_steps(r2bnd, right)
    else:
        r2_steps = _read_visible_start_steps(right, r2=True)
    r2_steps += _read_through_steps(r1bnd, left, "R2")
    plan = ReadPlan(r1bnd, r2bnd, r1_steps, r2_steps, r1_up, r2_up)
    return plan


def _read_visible_start_steps(arm, r2=False):
    """Read-start steps for a side WITHOUT a recognised boundary: the read begins
    at the outer end of its own arm, so every token on that arm is read-visible
    (R1 in written order; R2 bottom-strand in reversed, rc'd order)."""
    seq = reversed(arm) if r2 else arm
    return [_start_step(t, r2=r2) for t in seq]


# --- public entry point ---------------------------------------------------


def compile_read_plan(orientation, left, right, ctx, r1_primer=None, r2_primer=None):
    """Build ``(r1_mods, r2_mods, ReadPlan)`` or ``None`` (no primer found).

    Each side is resolved independently.  When only one side has a recognised
    primer the other side's modifiers are returned empty and flagged via
    ``plan.r1_is_plan``/``r2_is_plan`` so the caller fills the legacy per-token
    emission for the read-visible side.  Single-end with a recognised R1 site
    also starts downstream.
    """
    # Every side is resolved independently; a plan is built whenever ANY side
    # has a recognised primer, and BOTH reads are then routed through a coherent
    # arm-local plan (the unrecognised side uses all of its own read-visible
    # tokens plus the opposite correct binder's read-through).  A scheme with no
    # primer on either side (e.g. DBiT / TSO handle) resolves to ``None`` and
    # falls back to the legacy per-token emission.
    r1b = _resolve_boundary(left, "R1", r1_primer)
    r2b = _resolve_boundary(right, "R2", r2_primer)
    plan = _build_plan(left, right, r1b, r2b, orientation)
    if plan is None:
        return None
    mods1 = [_PlanExecutor(plan.r1_steps, 0, ctx.force_anywhere)]
    mods2 = [_PlanExecutor(plan.r2_steps, 1, ctx.force_anywhere)]
    plan.r1_own_arm_adapters = mods1[0].own_arm_adapters()
    plan.r2_own_arm_adapters = mods2[0].own_arm_adapters()
    return mods1, mods2, plan
