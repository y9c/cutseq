"""General ``-a`` scheme engine for cutseq (full-library grammar).

The legacy scheme is a fixed-slot string parsed by ``BarcodeConfig``. This
module is the (new) grammar engine: it uses the *same* element vocabulary
(adapters, inline barcodes, ``N`` captures, ``X`` masks) but takes a
space-separated, ordered full-library template (see grammar below).

Grammar (space-separated tokens, every character shell-safe)::

    AC.GT...     exact adapter / exact homopolymer   -> match & trim
    ac.gt...     inline barcode (lowercase)          -> match & capture
    NNNN  / N<n> random bases                        -> capture into read name
    XXXX  / X<n> mask                                -> trim (not captured)
    AAA...AAA    variable poly-A/T tail             -> trim whole 3' run (>= min)
    TTT...TTT
    : | + | -    insert (keep) + R1/R2 split        -> unstranded / sense / antisense

Everything left of the insert marker applies to R1 (in written order).
Everything right of it applies to R2 as the *reverse-complemented, reversed*
sequence. Capture labels come from a separate ``--barcode-names`` list (in
written order); if omitted they are numbered ``barcode1, barcode2, ...``.
"""

import re
from pathlib import Path

import cutadapt.modifiers as _mods
from cutadapt.adapters import (
    BackAdapter,
    PrefixAdapter,
    RightmostFrontAdapter,
    SuffixAdapter,
)

# --- parsing ---------------------------------------------------------------

# unanchored poly-tail form so a poly tail may be followed by further tokens
_POLY_RUN_RE = re.compile(r"([ACGT])\1*\.\.\.\1*")

MIN_POLY_LEN = 3  # minimum run of a homopolymer to be trimmed as a poly tail


class _Token:
    __slots__ = ("kind", "value", "label", "index", "options")

    def __init__(self, kind, value, label=None, options=None):
        self.kind = kind      # adp | inline | capture | mask | polytail
        self.value = value    # seq / length / (base,) as appropriate
        self.label = label    # capture label
        self.index = None     # written capture ordinal (fixes name tag order)
        self.options = options or {}  # per-part extras (max_errors, min_overlap, min_len...)


_INSERT = frozenset(":+-")


def _scan_polytail(s, i):
    """If a dot-delimited poly-A/T tail (``B...B``) starts at ``i``, return
    ``(end_index, base)``, else ``None``. Unanchored so a tail may be followed
    by further tokens (e.g. ``AAA...AAATTT``)."""
    m = _POLY_RUN_RE.match(s[i:])
    if m:
        return i + m.end(), m.group(1)
    return None


def _scan_run_or_num(s, i, ch):
    """Scan an ``N``/``X`` token from index ``i``.
    ``ch`` is one of ``N`` (capture) / ``X`` (mask). Returns a _Token.
    """
    j = i
    n = len(s)
    # numeric shorthand: N8 / X3 (digits follow the letter)
    if j + 1 < n and s[j + 1].isdigit() and not s[j + 1 : j + 2].isalpha():
        k = j + 1
        while k < n and s[k].isdigit():
            k += 1
        length = int(s[j + 1 : k])
        kind = "capture" if ch == "N" else "mask"
        return _Token(kind, length)
    while j < n and s[j] == ch:
        j += 1
    kind = "capture" if ch == "N" else "mask"
    return _Token(kind, j - i)


def tokenize(scheme):
    """Split a space-free scheme into a list of ``_Token`` objects.

    Tokens are delimited by character-class transitions (no whitespace
    required): uppercase runs -> adapter, lowercase runs -> inline barcode,
    ``N``/``X`` (+ optional length) -> capture/mask, ``: + -`` -> insert marker,
    and ``A...A``-style dot-delimited runs -> poly-A/T tail.
    """
    tokens = []
    i = 0
    n = len(scheme)
    while i < n:
        ch = scheme[i]
        if ch in _INSERT:
            tokens.append(_Token("insert", ch))
            i += 1
        elif ch in "Xx":
            tokens.append(_scan_run_or_num(scheme, i, "X"))
            i = _token_end(scheme, i, "X")
        elif ch in "Nn":
            tokens.append(_scan_run_or_num(scheme, i, "N"))
            i = _token_end(scheme, i, "N")
        elif ch in "ACGT":
            # Numeric shorthand for a homopolymer adapter: `A12` -> 12 As.
            if i + 1 < n and scheme[i + 1].isdigit():
                k = i + 1
                while k < n and scheme[k].isdigit():
                    k += 1
                tokens.append(_Token("adp", ch * int(scheme[i + 1:k])))
                i = k
                continue
            # A poly tail may start here directly (e.g. `AAA...AAAGGG`).
            pt = _scan_polytail(scheme, i)
            if pt is not None:
                end, base = pt
                tokens.append(_Token("polytail", base))
                i = end
                continue
            # Otherwise scan the maximal uppercase run.
            j = i
            while j < n and scheme[j] in "ACGT":
                j += 1
            # If the run is directly followed by a poly-tail (dot-delimited)
            # homopolymer, split it off: the trailing homopolymer run of the
            # preceding uppercase bases is the poly-tail's leading run (e.g.
            # `ACGTAAA...AAA` -> adapter `ACGT` + poly-A `AAA...AAA`).
            if j < n and scheme[j] == ".":
                base = scheme[j - 1]
                k = j
                while k > i and scheme[k - 1] == base:
                    k -= 1
                m = _POLY_RUN_RE.match(scheme[k:])
                if m:
                    if k > i:
                        tokens.append(_Token("adp", scheme[i:k]))
                    tokens.append(_Token("polytail", m.group(1)))
                    i = k + m.end()
                    continue
            tokens.append(_Token("adp", scheme[i:j]))
            i = j
        elif ch in "acgt":
            j = i
            while j < n and scheme[j] in "acgt":
                j += 1
            tokens.append(_Token("inline", scheme[i:j]))
            i = j
        elif ch == ">":
            # `>SEQUENCE` declares a 3' read-through (BackAdapter) adapter.
            j = i + 1
            while j < n and scheme[j] in "ACGT":
                j += 1
            tokens.append(_Token("back", scheme[i + 1:j]))
            i = j
        else:
            raise ValueError(f"cutseq: cannot parse scheme character {ch!r} at {i}")
    return tokens


def _token_end(scheme, i, ch):
    """Return the index just past the N/X token whose first char is at ``i``."""
    n = len(scheme)
    j = i + 1
    while j < n and (scheme[j].isdigit() or scheme[j] == ch):
        j += 1
    return j


def parse_scheme(scheme, auto_inline=True):
    """Parse a full-library scheme.

    Returns ``(orientation, left, right)``: *left* are tokens written before
    the insert marker (applied to R1), *right* are tokens written after it
    (applied to R2 as reverse-complemented, reversed sequence). *orientation*
    is one of ``':'``, ``'+'``, ``'-'`` or ``None`` (no insert marker).

    With ``auto_inline``, uppercase adapter runs that sit strictly between two
    known sequencing primers are reclassified as inline barcodes (a barcode the
    user wrote in uppercase by mistake).
    """
    orientation = None
    left, right = [], []
    side = left
    for tok in tokenize(scheme):
        if tok.kind == "insert":
            if orientation is not None:
                raise ValueError("cutseq: multiple insert markers in scheme")
            orientation = tok.value
            side = right
            continue
        side.append(tok)
    if auto_inline:
        _auto_detect_inline(left, right)
    return orientation, left, right


def _auto_detect_inline(left, right):
    """Reclassify/ split uppercase adapters into inline barcodes + primers.

    The outermost adapters of a library scheme are the sequencing primers (p5
    on the 5' side, p7 on the 3' side). A fixed uppercase sequence adjacent to
    a primer is an inline barcode the user wrote in uppercase by mistake — it
    merges with the primer into one run (e.g. ``ATCACGAGATCGGAAGAGCACACGTC``).
    This splits ``barcode + primer`` (right side) or ``primer + barcode``
    (left side) back into an ``inline`` token plus the primer, and reclassifies
    any standalone inner adapter as ``inline``. Only fires when the scheme's
    outermost adapters are recognized sequencing primers, so genuinely custom
    schemes are never altered.
    """
    import logging

    from .primers import (MIN_PRIMER_MATCH, SEQUENCING_PRIMERS, _norm, _rc)

    # A detected barcode must be at least this long; shorter leftovers are
    # partial-primer match noise, not real barcodes.
    MIN_INLINE_BC = 4

    def _all_adp():
        return [t.value for t in left + right if t.kind == "adp"]

    all_adp = _all_adp()
    if len(all_adp) < 2:
        return

    primers = list({_norm(p) for p in SEQUENCING_PRIMERS.values()})

    def _primer_overlap(v, strand):
        """Longest suffix of *strand* (>= MIN_PRIMER_MATCH) that is a prefix
        or suffix of *v*. Returns (primer_part, is_prefix_of_v) or None."""
        hi = min(len(v), len(strand))
        for L in range(hi, MIN_PRIMER_MATCH - 1, -1):
            if v[:L] == strand[-L:]:
                return strand[-L:], True
            if v[-L:] == strand[-L:]:
                return strand[-L:], False
        return None

    def _is_pure_primer(value):
        """True if the whole value is a terminal fragment of a known primer
        (>= MIN_PRIMER_MATCH bp), so it is a primer itself, not a barcode."""
        v = _norm(value)
        for p in primers:
            for strand in (p, _rc(p)):
                if len(strand) < MIN_PRIMER_MATCH:
                    continue
                if len(v) >= MIN_PRIMER_MATCH and \
                   (strand.endswith(v) or strand.startswith(v)):
                    return True
        return False

    def _match_terminal(value):
        """Return (primer, barcode) if *value* = primer+barcode (primer at the
        value's 5' end) or barcode+primer (primer at the 3' end), with the
        barcode on the other side being non-trivial; else None.

        The longest primer match across all known primers wins, and a value
        that is itself a known primer is never split (its barcode side is
        empty).
        """
        v = _norm(value)
        if _is_pure_primer(v):
            return None
        best = None
        for p in primers:
            for strand in (p, _rc(p)):
                hit = _primer_overlap(v, strand)
                if hit is None:
                    continue
                primer_part, at_5prime = hit
                bc = v[len(primer_part):] if at_5prime else v[:len(v) - len(primer_part)]
                if not (MIN_INLINE_BC <= len(bc) < len(primer_part)):
                    continue
                if best is None or len(primer_part) > len(best[0]):
                    best = (primer_part, bc)
        return best

    def _as_inline(tok):
        return _Token("inline", tok.value, tok.label, options=tok.options)

    def _fix_side(side, outer_is_first):
        adp_idx = [i for i, t in enumerate(side) if t.kind == "adp"]
        if not adp_idx:
            return
        outer_i = adp_idx[0] if outer_is_first else adp_idx[-1]
        outer = side[outer_i]
        # Only touch a side whose outermost adapter is a standard sequencing
        # primer (pure or merged with a barcode); custom schemes are left as-is.
        if not (_is_pure_primer(outer.value) or _match_terminal(outer.value)):
            return
        # Split a merged barcode off the outermost adapter first.
        hit = _match_terminal(outer.value)
        if hit is not None:
            primer, bc = hit
            logging.warning(
                "Auto-detected inline barcode %r adjacent to sequencing "
                "primer %r; treating it as an inline barcode.",
                bc, primer,
            )
            outer.value = primer
            outer.kind = "adp"
            # Left side: primer+barcode -> barcode after primer.
            # Right side: barcode+primer -> barcode before primer.
            insert_at = outer_i + (1 if outer_is_first else 0)
            side.insert(insert_at,
                        _Token("inline", bc, outer.label, options=outer.options))
            adp_idx = [i for i, t in enumerate(side) if t.kind == "adp"]
        # Reclassify remaining standalone inner adapters as inline barcodes.
        if len(adp_idx) < 2:
            return
        for i in adp_idx:
            if (outer_is_first and i == adp_idx[0]) or \
               (not outer_is_first and i == adp_idx[-1]):
                continue
            logging.warning(
                "Auto-detected inline barcode %r between sequencing primers; "
                "treating it as an inline barcode (write it lowercase if "
                "this is unintended).",
                side[i].value,
            )
            side[i] = _as_inline(side[i])

    _fix_side(left, outer_is_first=True)
    _fix_side(right, outer_is_first=False)


def assign_labels(tokens, names):
    """Assign a capture label and written ordinal to each capture/inline token."""
    names = list(names) if names else []
    n = 0
    for t in tokens:
        if t.kind in ("capture", "inline"):
            t.index = n
            t.label = names[n] if n < len(names) else f"barcode{n + 1}"
            n += 1


# --- YAML-part scheme (``-A file.yaml``) -----------------------------------

_OPT_KEYS = {
    "max_errors", "min_overlap", "min_len",
    "read_wildcards", "adapter_wildcards", "indels", "force_anywhere",
}


def parse_scheme_parts(data):
    """Parse a YAML scheme into ``(orientation, left, right)``.

    ``data`` is the decoded document; it may be a bare list of parts, or a dict
    with a ``parts`` list. Parts look like::

        - type: adapter
          seq: AGTTCTACAGTCCGACGATC
          max_errors: 0.2
        - type: insert        # R1/R2 split marker: value is '+', '-' or ':'
          value: "+"
        - type: adapter
          seq: AGATCGGAAGAGCACACGTC

    ``type`` maps to the grammar kinds: ``adapter`` -> adp, ``inline`` -> inline,
    ``capture`` -> capture (N), ``mask`` -> mask (X), ``polytail`` -> polytail,
    ``insert`` -> R1/R2 insert marker. For ``adapter``/``inline`` a ``seq`` field
    is required; for ``capture``/``mask`` a ``length``; for ``polytail`` a
    ``base``; for ``insert`` a ``value`` (one of ``+``/``-``/``:``). Per-part
    settings -- ``max_errors``, ``min_overlap``, ``min_len``, ``read_wildcards``,
    ``adapter_wildcards``, ``indels``, ``force_anywhere`` -- are carried through
    to the modifier.

    Parts before the ``insert`` part belong to R1 (left); parts after it belong
    to R2 (right). Without an ``insert`` part, all parts apply to a single read.
    """
    if isinstance(data, dict):
        raw_parts = data.get("parts")
    else:
        raw_parts = data
    if not isinstance(raw_parts, list) or not raw_parts:
        raise ValueError("cutseq: YAML scheme needs a non-empty 'parts' list")

    orientation = None
    left, right = [], []
    side = left
    for i, part in enumerate(raw_parts):
        typ = part.get("type")
        if typ is None:
            raise ValueError(f"cutseq: part #{i} is missing 'type'")
        if typ == "insert":
            if orientation is not None:
                raise ValueError(f"cutseq: multiple 'insert' parts (at #{i})")
            orientation = part.get("value")
            if orientation not in _INSERT:
                raise ValueError(
                    f"cutseq: insert part #{i} 'value' must be '+', '-' or ':'"
                )
            side = right
            continue
        options = {k: part[k] for k in _OPT_KEYS if k in part}
        label = part.get("label")
        if typ == "adapter":
            seq = part.get("seq")
            if not seq:
                raise ValueError(f"cutseq: adapter part #{i} needs 'seq'")
            tok = _Token("adp", seq, label, options)
        elif typ == "adapter_back":
            seq = part.get("seq")
            if not seq:
                raise ValueError(f"cutseq: adapter_back part #{i} needs 'seq'")
            tok = _Token("back", seq, label, options)
        elif typ == "inline":
            seq = part.get("seq")
            if not seq:
                raise ValueError(f"cutseq: inline part #{i} needs 'seq'")
            tok = _Token("inline", seq, label, options)
        elif typ == "capture":
            length = part.get("length")
            if not length or length < 1:
                raise ValueError(f"cutseq: capture part #{i} needs positive 'length'")
            tok = _Token("capture", int(length), label, options)
        elif typ == "mask":
            length = part.get("length")
            if not length or length < 1:
                raise ValueError(f"cutseq: mask part #{i} needs positive 'length'")
            tok = _Token("mask", int(length), label, options)
        elif typ == "polytail":
            base = part.get("base")
            if not base:
                raise ValueError(f"cutseq: polytail part #{i} needs 'base'")
            tok = _Token("polytail", base, None, options)
        else:
            raise ValueError(f"cutseq: unknown part type {typ!r} at #{i}")
        side.append(tok)
    return orientation, left, right


# --- compile: lazy token graph -> cutadapt-native modifiers -----------------


_ADAPTER_KW = ("max_errors", "min_overlap", "read_wildcards",
               "adapter_wildcards", "indels", "force_anywhere")


class _ConditionalCutter(_mods.SingleEndModifier):
    """Like legacy ``ConditionalCutter``: trim n bases only if an adapter was
    matched or the read is long enough. cutadapt 5.2 has no native equivalent."""

    def __init__(self, length, force_trim_min_length=50, conditional=True,
                 capture=False):
        self.length = length
        self.force_trim_min_length = force_trim_min_length
        self.conditional = conditional
        self.capture = capture

    def __call__(self, read, info):
        if self.conditional and not info.matches and len(read.sequence) < self.force_trim_min_length:
            return read
        if self.length > 0:
            if self.capture:
                info.cut_prefix = read.sequence[: self.length]
            return read[self.length:]
        if self.capture:
            info.cut_suffix = read.sequence[self.length:]
        return read[: self.length]

    def __repr__(self):
        return f"ConditionalCutter(length={self.length})"


def _adapter_kw(options):
    return {k: options[k] for k in _ADAPTER_KW
            if k in options and options[k] is not None}


class _Ctx:
    """Compilation context shared by token emitters."""

    __slots__ = ("conditional", "force_trim_min_length", "force_anywhere")

    def __init__(self, conditional=True, force_trim_min_length=50, force_anywhere=False):
        self.conditional = conditional
        self.force_trim_min_length = force_trim_min_length
        self.force_anywhere = force_anywhere


def _emit_adp_five(tok, ctx):
    kw = {"max_errors": 0.2, "min_overlap": 10, **_adapter_kw(tok.options)}
    return _mods.AdapterCutter([RightmostFrontAdapter(tok.value, **kw)])


def _emit_adp_three(tok, ctx):
    kw = {"max_errors": 0.2, "min_overlap": 3,
          "force_anywhere": ctx.force_anywhere, **_adapter_kw(tok.options)}
    return _mods.AdapterCutter([BackAdapter(tok.value, **kw)])


def _emit_back_five(tok, ctx):
    return _mods.AdapterCutter([BackAdapter(tok.value, **_adapter_kw(tok.options))])


def _emit_inline_five(tok, ctx):
    kw = {"max_errors": 0.2, **_adapter_kw(tok.options)}
    mod = _mods.AdapterCutter([PrefixAdapter(tok.value, **kw)])
    mod._is_inline = True
    return mod


def _emit_inline_three(tok, ctx):
    return _ConditionalCutter(-len(tok.value), ctx.force_trim_min_length,
                              conditional=False)


def _emit_inline_se_three(tok, ctx):
    mod = _mods.AdapterCutter([SuffixAdapter(tok.value, max_errors=0.2)])
    mod._is_inline = True
    return mod


def _emit_capture_five(tok, ctx):
    return _ConditionalCutter(tok.value, ctx.force_trim_min_length,
                              conditional=False, capture=True)


def _emit_capture_three(tok, ctx):
    return _ConditionalCutter(-tok.value, ctx.force_trim_min_length,
                              conditional=ctx.conditional, capture=True)


def _emit_capture_se_three(tok, ctx):
    return _ConditionalCutter(-tok.value, ctx.force_trim_min_length,
                              conditional=False, capture=True)


def _emit_mask_five(tok, ctx):
    return _ConditionalCutter(tok.value, ctx.force_trim_min_length,
                              conditional=False)


def _emit_mask_three(tok, ctx):
    return _ConditionalCutter(-tok.value, ctx.force_trim_min_length,
                              conditional=ctx.conditional)


def _emit_mask_se_three(tok, ctx):
    return _ConditionalCutter(-tok.value, ctx.force_trim_min_length,
                              conditional=False)


def _emit_polytail_five(tok, ctx):
    return PolyTailModifier(tok.value, tok.options.get("min_len") or MIN_POLY_LEN)


# Token-kind registry: kind -> (five_emitter, three_emitter, se_three_emitter).
#   five       = written-side 5' emission (FrontAdapter / positive cuts)
#   three      = paired-end mirrored 3' emission (BackAdapter / negative cuts)
#   se_three   = single-end 3' emission on the same read
# A None emitter means the token is not emitted on that side. Adding a new
# token kind is one registry entry plus one phase entry below.
_EMITTERS = {
    "adp": (_emit_adp_five, _emit_adp_three, _emit_adp_three),
    "back": (_emit_back_five, None, None),
    "inline": (_emit_inline_five, _emit_inline_three, _emit_inline_se_three),
    "capture": (_emit_capture_five, _emit_capture_three, _emit_capture_se_three),
    "mask": (_emit_mask_five, _emit_mask_three, _emit_mask_se_three),
    "polytail": (_emit_polytail_five, _emit_polytail_five, _emit_polytail_five),
}

# Paired-end phase order (legacy pipeline order). "by_end" pairs adapters on
# both reads (5' front then 3' read-through); "own" emits the read's own 3'
# read-through; "side" emits written-side (left) tokens on both reads before
# right-side tokens, inside each inline/capture/mask phase.
_PAIRED_PHASES = (
    ("adp", "by_end"),
    ("back", "own"),
    ("inline", "side"),
    ("capture", "side"),
    ("mask", "side"),
    ("polytail", "side"),
)

_SINGLE_KINDS = ("adp", "back", "inline", "capture", "mask", "polytail")


def _side(tokens, kind):
    return [t for t in tokens if t.kind == kind]


def compile_tokens(orientation, left, right, paired=True, conditional_cutter=True,
                   force_trim_min_length=50, force_anywhere=False):
    """Compile the parsed token graph into ``(r1_mods, r2_mods)``.

    Each token describes one physical end of the molecule, so it applies to
    both reads: written-side at the 5' end (FrontAdapter / positive
    unconditional cuts) and mirrored at the other read's 3' end (read-through
    BackAdapter / negative conditional cuts). Phase order matches legacy:
    adapters (front then read-through), explicit back, inline, captures, then
    masks — with written-side (left) tokens before mirrored (right) tokens
    inside the inline/capture/mask phases.
    """
    ctx = _Ctx(conditional_cutter, force_trim_min_length, force_anywhere)
    rev_left, rev_right = list(reversed(left)), list(reversed(right))

    if not paired or orientation is None:
        r1 = []
        for kind in _SINGLE_KINDS:
            five, _, se_three = _EMITTERS[kind]
            r1 += [five(t, ctx) for t in _side(left, kind) if five is not None]
            r1 += [se_three(t, ctx) for t in _side(rev_right, kind)
                   if se_three is not None]
        return r1, []

    mods1, mods2 = [], []
    for kind, shape in _PAIRED_PHASES:
        five, three, _ = _EMITTERS[kind]
        if shape == "by_end":
            if five:
                mods1 += [five(t, ctx) for t in _side(left, kind)]
                mods2 += [five(_for_r2(t), ctx) for t in _side(rev_right, kind)]
            if three:
                mods1 += [three(t, ctx) for t in _side(rev_right, kind)]
                mods2 += [three(_for_r2(t), ctx) for t in _side(rev_left, kind)]
        elif shape == "own":
            if five:
                mods1 += [five(t, ctx) for t in _side(left, kind)]
                mods2 += [five(_for_r2(t), ctx) for t in _side(rev_right, kind)]
        else:
            if five:
                mods1 += [five(t, ctx) for t in _side(left, kind)]
            if three:
                mods2 += [three(_for_r2(t), ctx) for t in _side(rev_left, kind)]
                mods1 += [three(t, ctx) for t in _side(rev_right, kind)]
            if five:
                mods2 += [five(_for_r2(t), ctx) for t in _side(rev_right, kind)]
    return mods1, mods2


def _se_id_name(read, info):
    """``{id}``: strip any trailing comment (identical to cutadapt Renamer)."""
    read.name = read.name.split(maxsplit=1)[0]
    return read


def _make_fast_renamer(paired, has_captures, name_format, capture_separator):
    """Build a lean renamer that computes only the variables its template needs.

    cutadapt's generic ``Renamer``/``PairedEndRenamer`` builds a full variable
    dict (comment, header, adapter_name, match_sequence, ...) plus a
    ``SimpleNamespace`` per read on every call even when the template uses just
    ``id``/``cut_prefix``/``cut_suffix``. The default cutseq templates only use
    those, so we special-case them and fall back to the native renamer for any
    custom ``--name-format`` template. Output is byte-identical to the native
    renamer for the same template.
    """
    if name_format is not None:
        return None
    if not has_captures:
        tpl = "{id}"
        if paired:
            def _rename(read1, read2, info1, info2):
                read1.name = read1.name.split(maxsplit=1)[0]
                read2.name = read2.name.split(maxsplit=1)[0]
                return read1, read2
            _rename._template = tpl
            return _rename
        _se_id_name._template = tpl
        return _se_id_name
    if paired:
        if capture_separator is not None:
            sep = capture_separator
            tpl = f"{{id}}{sep}{{r1.cut_prefix}}{sep}{{r2.cut_prefix}}"

            def _rename(read1, read2, info1, info2):
                id1 = read1.name.split(maxsplit=1)[0]
                id2 = read2.name.split(maxsplit=1)[0]
                p1 = info1.cut_prefix if info1.cut_prefix else ""
                p2 = info2.cut_prefix if info2.cut_prefix else ""
                read1.name = f"{id1}{sep}{p1}{sep}{p2}"
                read2.name = f"{id2}{sep}{p1}{sep}{p2}"
                return read1, read2
        else:
            tpl = "{id}_{r1.cut_prefix}{r2.cut_prefix}"

            def _rename(read1, read2, info1, info2):
                id1 = read1.name.split(maxsplit=1)[0]
                id2 = read2.name.split(maxsplit=1)[0]
                p1 = info1.cut_prefix if info1.cut_prefix else ""
                p2 = info2.cut_prefix if info2.cut_prefix else ""
                read1.name = f"{id1}_{p1}{p2}"
                read2.name = f"{id2}_{p1}{p2}"
                return read1, read2
        _rename._template = tpl
        return _rename
    if capture_separator is not None:
        sep = capture_separator
        tpl = f"{{id}}{sep}{{cut_prefix}}{sep}{{cut_suffix}}"

        def _rename(read, info):
            p = info.cut_prefix if info.cut_prefix else ""
            s = info.cut_suffix if info.cut_suffix else ""
            read.name = f"{read.name.split(maxsplit=1)[0]}{sep}{p}{sep}{s}"
            return read
    else:
        tpl = "{id}_{cut_prefix}{cut_suffix}"

        def _rename(read, info):
            p = info.cut_prefix if info.cut_prefix else ""
            s = info.cut_suffix if info.cut_suffix else ""
            read.name = f"{read.name.split(maxsplit=1)[0]}_{p}{s}"
            return read
    _rename._template = tpl
    return _rename


class _FastPairedEndRenamer(_mods.PairedEndModifier):
    """Wraps a fast paired rename callable so it plugs into the pipeline."""

    def __init__(self, fn):
        self._fn = fn
        self._template = fn._template

    def __call__(self, read1, read2, info1, info2):
        return self._fn(read1, read2, info1, info2)


def make_renamer(paired, has_captures=False, name_format=None,
                 capture_separator=None):
    """Native cutadapt renamer. Appends captured ``N`` UMIs to the read name
    only when the scheme declares capture tokens (legacy naming). Single-end
    UMIs sit at the 3' end, so both cut_prefix and cut_suffix are appended.

    ``name_format`` overrides the template entirely (cutadapt brace syntax,
    e.g. ``{id}_{cut_prefix}{cut_suffix}``). ``capture_separator`` inserts a
    separator between captures in the default template. Defaults reproduce
    legacy naming exactly.
    """
    if name_format is not None:
        tpl = name_format
    elif not has_captures:
        tpl = "{id}"
    elif paired:
        if capture_separator is not None:
            tpl = ("{id}" + capture_separator + "{r1.cut_prefix}"
                   + capture_separator + "{r2.cut_prefix}")
        else:
            tpl = "{id}_{r1.cut_prefix}{r2.cut_prefix}"
    elif capture_separator is not None:
        tpl = ("{id}" + capture_separator + "{cut_prefix}"
               + capture_separator + "{cut_suffix}")
    else:
        tpl = "{id}_{cut_prefix}{cut_suffix}"
    fast = _make_fast_renamer(paired, has_captures, name_format,
                              capture_separator)
    if fast is not None:
        if paired:
            return _FastPairedEndRenamer(fast)
        return fast
    return _mods.PairedEndRenamer(tpl) if paired else _mods.Renamer(tpl)


class CompiledScheme:
    """A parsed library scheme that compiles lazily to cutadapt modifiers.

    Holds the token graph (``orientation``, ``left``, ``right``) plus the
    compile-affecting settings. ``modifiers()`` builds the cutadapt modifier
    lists for the single/paired pipeline on demand; ``renamer()`` builds the
    native name renamer. Adding a token kind is one entry in ``_EMITTERS``
    plus one phase entry — no changes needed here.
    """

    def __init__(self, orientation, left, right, conditional_cutter=True,
                 force_trim_min_length=50, force_anywhere=False):
        self.orientation = orientation
        self.left = left
        self.right = right
        self.conditional_cutter = conditional_cutter
        self.force_trim_min_length = force_trim_min_length
        self.force_anywhere = force_anywhere

    @property
    def has_captures(self):
        return any(t.kind == "capture" for t in self.left + self.right)

    @property
    def has_inline(self):
        return any(t.kind == "inline" for t in self.left + self.right)

    def modifiers(self, paired=True):
        """Return ``(r1_mods, r2_mods)`` lists of cutadapt modifier callables."""
        cache = getattr(self, "_mods_cache", {})
        if paired in cache:
            return cache[paired]
        mods1, mods2 = compile_tokens(
            self.orientation, self.left, self.right, paired=paired,
            conditional_cutter=self.conditional_cutter,
            force_trim_min_length=self.force_trim_min_length,
            force_anywhere=self.force_anywhere,
        )
        # Cache the inline-barcode Adapter objects from THIS compile so a
        # subsequent inline_adapters() call returns the exact same instances
        # the pipeline's AdapterCutter steps use (Adapter equality is
        # identity-based, so recompiling would break --ensure-inline-barcode).
        cache[paired] = (mods1, mods2)
        self._inline_cache = getattr(self, "_inline_cache", {})
        self._inline_cache[paired] = self._collect_inline(mods1, mods2)
        return mods1, mods2

    @staticmethod
    def _collect_inline(mods1, mods2):
        def _grab(mods):
            out = []
            for m in mods:
                if getattr(m, "_is_inline", False):
                    out.extend(getattr(m, "adapters", ()))
            return out

        return _grab(mods1), _grab(mods2)

    def renamer(self, paired=True, name_format=None, capture_separator=None):
        """Return a native cutadapt renamer (see ``make_renamer``)."""
        return make_renamer(paired, self.has_captures, name_format,
                            capture_separator)

    def inline_adapters(self, paired=True):
        """Return ``(r1_adps, r2_adps)`` — the inline-barcode Adapter objects
        used by the compiled 5' AdapterCutter steps (empty lists if the scheme
        has no inline barcodes). These are the exact objects recorded in
        ``info.matches``, so they can drive ``--ensure-inline-barcode``
        filtering without re-deriving sequence identity."""
        cache = getattr(self, "_inline_cache", {})
        if paired not in cache:
            self.modifiers(paired)
        return cache[paired]


def rc(seq):
    return _rc(seq)


def _for_r2(tok):
    """Return the token as it should be applied to R2 (reverse-complemented).

    R2 is the RC of the top-strand 3' side, so adapters/inline barcodes are
    reverse-complemented. Captures (N) and masks (X) are length-based and RC is
    the identity for them; poly tails are unchanged.
    """
    if tok.kind in ("adp", "back", "inline"):
        return _Token(tok.kind, _rc(tok.value), tok.label, options=tok.options)
    return tok


def _rc(seq):
    return str(seq).translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


class PolyTailModifier(_mods.SingleEndModifier):
    """Trim a 3' homopolymer run of *base* if it is at least *min_len* long."""

    def __init__(self, base, min_len=MIN_POLY_LEN):
        self.base = base
        self.min_len = min_len

    def __call__(self, read, info):
        seq = read.sequence
        i = len(seq)
        while i > 0 and seq[i - 1] == self.base:
            i -= 1
        if len(seq) - i >= self.min_len:
            return read[:i]
        return read

    def __repr__(self):
        return f"PolyTail({self.base},min={self.min_len})"


def build_modifiers(scheme, barcode_names=None, paired=True,
                    conditional_cutter=True, force_trim_min_length=50):
    """Compile a library-grammar scheme into ``(r1_mods, r2_mods, orientation)``."""
    orientation, left, right = parse_scheme(scheme)
    return build_modifiers_from_parts(
        orientation, left, right, barcode_names, paired,
        conditional_cutter, force_trim_min_length,
    )


def build_modifiers_from_parts(orientation, left, right, barcode_names=None,
                               paired=True, conditional_cutter=True,
                               force_trim_min_length=50):
    """Compile parsed part tokens into ``(r1_mods, r2_mods, orientation)``."""
    assign_labels(left + right, barcode_names)
    r1, r2 = compile_tokens(orientation, left, right, paired,
                            conditional_cutter, force_trim_min_length)
    return r1, r2, orientation


def load_scheme_file(path):
    """Load a ``-A`` scheme value: a built-in / inline string or a YAML file.

    Returns ``(is_yaml, value)`` where ``value`` is either the scheme string or
    the ``(orientation, left, right)`` token triple.
    """
    p = path if isinstance(path, (str, bytes)) else ""
    low = str(p).lower()
    if low.endswith((".yaml", ".yml")) or (isinstance(p, str) and Path(p).exists()):
        try:
            import yaml
        except ImportError:  # pragma: no cover
            raise ValueError("cutseq: YAML scheme files require 'pyyaml'")
        with open(p, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
        return True, parse_scheme_parts(data)
    return False, path
