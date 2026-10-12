#!/usr/bin/env python3
"""Generate docs/adapters.html as a self-contained interactive adapter viewer.

The emitted page is a real ``.html`` file (no Jekyll front-matter, full ``<html>``
shell, inline <script> + <style>, inline ``adapdata`` JSON). Jekyll / GitHub Pages
copies ``.html`` verbatim and never runs it through Kramdown / jekyll-spaceship, so
the embedded JavaScript is never corrupted.

Every built-in scheme from cutseq/adapters.toml is tokenized with cutseq's own
grammar (cutseq.grammar.tokenize), so each segment is colored and labelled by
its REAL role (adapter / UMI capture / mask / inline barcode / poly-tail /
read-through / insert marker). The page embeds all schemes as JSON and renders
an interactive, searchable, filterable gallery:
  * legend of token roles and their colors
  * free-text search (name / description / scheme)
  * per-segment hover tooltip (role, length, capture index)
  * click a segment to copy just that segment
  * copy the whole scheme and the `cutseq -A NAME` command per scheme
"""
import html
import json
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

try:
    from cutseq.grammar import tokenize
except ImportError as e:
    sys.exit(f"Error importing tokenize from cutseq.grammar: {e}")

try:
    import tomllib  # Python 3.11+
except ImportError:
    try:
        import tomli as tomllib  # Python < 3.11
    except ImportError:
        sys.exit("Error: 'tomli' required for Python < 3.11.")

# --- token role metadata: kind -> (label, css class) -------------------------
ROLES = {
    "adp":      ("adapter / primer", "adp"),
    "inline":   ("inline barcode", "inline"),
    "capture":  ("UMI / capture", "capture"),
    "mask":     ("mask (trimmed, not kept)", "mask"),
    "polytail": ("homopolymer tail", "polytail"),
    "back":     ("3' read-through adapter", "back"),
    "insert":   ("insert marker (strand)", "insert"),
}

_INSERT = {"+": "sense", "-": "antisense", ":": "unstranded"}


def _display_len(kind, value):
    """Display run for captures/masks, capped so rows stay readable."""
    ch = {"capture": "N", "mask": "X"}[kind]
    return ch * min(value, 96)


# --- name each scheme token by its real role / known primer ----------------
try:
    from cutseq.primers import SEQUENCING_PRIMERS as _SP, _norm as _norm_seq
    _PRI_EXACT = {_norm_seq(v): k for k, v in _SP.items()}
    _P5 = _SP.get("P5 (flowcell)", "AATGATACGGCGACCACCGAGATCTACAC")
    _P7 = _SP.get("P7 (flowcell)", "CAAGCAGAAGACGGCATACGAGAT")
except Exception:
    _PRI_EXACT = {}
    _P5 = "AATGATACGGCGACCACCGAGATCTACAC"
    _P7 = "CAAGCAGAAGACGGCATACGAGAT"


def _tok_name(tok, side):
    """A human name for a scheme token (role-aware, known-primer aware)."""
    kind = tok.kind
    if kind == "capture":
        return f"UMI/capture N{int(tok.value)}" + (f" (rename {{{tok.index}}})" if tok.index else "")
    if kind == "mask":
        return f"mask X{int(tok.value)}"
    if kind == "insert":
        return f"insert marker ({_INSERT.get(tok.value, tok.value)})"
    if kind == "inline":
        return "inline barcode"
    if kind == "polytail":
        return f"poly-{tok.value} tail"
    if kind == "back":
        return "3' read-through adapter"
    # adapter / primer: match a known sequencing primer (exact, then terminal fragment)
    if tok.kind == "adp":
        if _PRI_EXACT:
            name = _PRI_EXACT.get(_norm_seq(tok.value))
            if name:
                return name
        try:
            from cutseq.primers import _terminal_match
            for k, v in sorted(_SP.items(), key=lambda kv: -len(kv[1])):
                if len(tok.value) >= 15 and _terminal_match(tok.value, v):
                    return k
        except Exception:
            pass
    return "5' primer/adapter (R1)" if side == "R1" else "3' primer/adapter (R2)"


def _token_json(tok, ordinal, name):
    """A token -> {k, s, n, nm} for the viewer (k=kind, s=display, n=capture ord, nm=name)."""
    kind = tok.kind
    base = {"k": kind, "nm": name}
    if kind == "capture":
        base.update({"s": _display_len(kind, tok.value), "n": ordinal})
    elif kind == "mask":
        base.update({"s": _display_len(kind, tok.value), "n": None})
    elif kind == "insert":
        base.update({"s": tok.value, "n": None})
    elif kind == "polytail":
        base.update({"s": f"{tok.value}...{tok.value}", "n": None})
    else:
        base.update({"s": tok.value, "n": None})
    return base


INSERT_PH = 30     # nominal length of the (unknown) library insert placeholder
POLY_RUN = 8       # nominal run length shown for homopolymer tails


def _scheme_json(name, info, tokens, construct):
    """One scheme -> the payload the viewer renderer consumes."""
    return {
        "name": name,
        "desc": info.get("description_name", ""),
        "scheme": info["scheme"],
        "points": info.get("points", []),
        "alias": info.get("alias_of"),
        "tokens": tokens,
        "construct": construct,
    }


def _build_construct(tokens, names):
    """Assemble the FULL top-strand molecule 5'->3'.

    = [P5 flowcell][i5 index] + R1-side + insert + R2-side + [i7 index][P7 flowcell]
    Every token becomes an annotated feature with 1-based coords; the flowcell
    anchors and index placeholders are added so the construct is full length.
    """
    side = "R1"
    specs = []  # (name, text, role, strand, side)
    specs.append(("P5 (flowcell)", _P5, "adp", "+", "R1"))
    specs.append(("i5 index (8N)", "N" * 8, "mask", "+", "R1"))
    for i, tok in enumerate(tokens):
        kind = tok.kind
        nm = names[i] if i < len(names) else None
        if kind == "insert":
            specs.append(("library insert (placeholder)", "N" * INSERT_PH, "insert", "+", side))
            side = "R2"
            continue
        if kind == "capture":
            specs.append((nm or f"UMI/capture N{int(tok.value)}", "N" * int(tok.value), "capture", "+", side))
        elif kind == "mask":
            specs.append((nm or f"mask X{int(tok.value)}", "N" * int(tok.value), "mask", "+", side))
        elif kind == "adp":
            specs.append((nm or ("5' primer/adapter (R1)" if side == "R1" else "3' primer/adapter (R2)"), tok.value, "adp", "+", side))
        elif kind == "inline":
            specs.append((nm or "inline barcode", tok.value, "inline", "+", side))
        elif kind == "polytail":
            specs.append((nm or f"poly-{tok.value} tail", tok.value * POLY_RUN, "polytail", "+", side))
        elif kind == "back":
            specs.append((nm or "3' read-through adapter", tok.value, "back", "-", side))
    specs.append(("i7 index (8N)", "N" * 8, "mask", "+", "R2"))
    specs.append(("P7 (flowcell)", _P7, "adp", "+", "R2"))

    parts, feats, cursor = [], [], 0
    for nm, text, role, strand, sd in specs:
        start, end = cursor + 1, cursor + len(text)
        feats.append({"name": nm, "type": "misc_feature", "role": role, "strand": strand,
                      "side": sd, "start": start, "end": end, "len": len(text), "seq": text})
        parts.append(text)
        cursor += len(text)
    return {"seq": "".join(parts), "size": cursor, "features": feats}


def _load_schemes():
    data = tomllib.loads(Path(project_root / "cutseq" / "adapters.toml").read_text("utf-8"))
    schemes = []
    for name, info in data.items():
        if not isinstance(info, dict) or "scheme" not in info:
            continue
        if not info.get("description_name"):
            continue
        # captures (UMI N) and inline barcodes are counted together for --rename {N}
        ordinal = 0
        raw_tokens = tokenize(info["scheme"])
        tokens = []
        names = []
        side = "R1"
        for tok in raw_tokens:
            if tok.kind in ("capture", "inline"):
                ordinal += 1
                tok.index = ordinal
            names.append(_tok_name(tok, side))
            if tok.kind == "insert":
                side = "R2"
        for i, tok in enumerate(raw_tokens):
            tok_ord = tok.index if tok.kind in ("capture", "inline") else None
            tokens.append(_token_json(tok, tok_ord, names[i]))
        construct = _build_construct(raw_tokens, names)
        schemes.append(_scheme_json(name, info, tokens, construct))
    return schemes


def _legend_html():
    rows = []
    for kind in ("adp", "capture", "mask", "inline", "polytail", "back", "insert"):
        label, css = ROLES[kind]
        rows.append(
            '<div class="adleg"><span class="leg-swatch seg-' + css + '"></span>'
            '<span class="leg-label">' + html.escape(label) + "</span></div>"
        )
    return '<div class="ad-legend">' + "".join(rows) + "</div>"


# --- viewer CSS + JS (kept free of Jekyll/Liquid `{{` sequences) -------------
_CSS = """
/* The viewer couples to the shared shell tokens (--cs-*) so light AND dark mode
   stay coherent with the rest of the site; the hardcoded fallbacks keep it usable
   even if a token is missing. */
.ad-toolbar{display:flex;gap:.75rem;align-items:center;flex-wrap:wrap;margin:1.1rem 0 .4rem;padding:.6rem .7rem;background:var(--cs-bg-soft,#f8fafc);border:1px solid var(--cs-border,#e2e8f0);border-radius:var(--cs-radius,10px);box-shadow:var(--cs-shadow-sm,0 1px 2px rgba(15,23,42,.06))}
.ad-search{flex:1 1 240px;font-size:.92rem;padding:.5rem .7rem;border:1px solid var(--cs-border-strong,#cbd5e1);border-radius:var(--cs-radius-sm,8px);min-width:190px;background:var(--cs-bg,#fff);color:var(--cs-text,#1f2937);outline:none;transition:border-color .12s,box-shadow .12s}
.ad-search::placeholder{color:var(--cs-faint,#94a3b8)}
.ad-search:focus{outline:none;border-color:var(--cs-primary,#0d9488);box-shadow:0 0 0 3px var(--cs-primary-ghost-2,rgba(13,148,136,.14))}
.ad-count{font-size:.82rem;color:var(--cs-muted,#64748b);white-space:nowrap;font-variant-numeric:tabular-nums}
.ad-legend{display:flex;gap:1rem;flex-wrap:wrap;align-items:center;margin:.1rem 0 1.1rem;padding:.6rem 1rem;background:var(--cs-bg-soft,#f8fafc);border:1px solid var(--cs-border,#e2e8f0);border-radius:var(--cs-radius,10px)}
.adleg{display:flex;align-items:center;gap:.42rem;font-size:.8rem;color:var(--cs-muted,#64748b);font-weight:500}
.leg-swatch{width:13px;height:13px;border-radius:4px;display:inline-block}
.adapter-card{margin:.7rem 0;border:1px solid var(--cs-border,#e2e8f0);border-radius:var(--cs-radius,10px);background:var(--cs-bg,#fff);overflow:hidden;box-shadow:var(--cs-shadow-sm,0 1px 2px rgba(15,23,42,.04));transition:border-color .15s,box-shadow .15s}
.adapter-card[open]{border-color:var(--cs-border-strong,#cbd5e1);box-shadow:var(--cs-shadow-md,0 6px 20px rgba(15,23,42,.08))}
.adapter-card>summary{display:flex;gap:.7rem;align-items:baseline;cursor:pointer;padding:.6rem .85rem;list-style:none;flex-wrap:wrap;background:var(--cs-bg-soft,#f8fafc)}
.adapter-card>summary:hover{background:var(--cs-primary-ghost,rgba(13,148,136,.05))}
.adapter-card>summary::-webkit-details-marker{display:none}
.adapter-card[data-alias]{border-left:4px solid var(--cs-accent-soft,#f59e0b)}
.card-name{font-family:var(--cs-font-mono,ui-monospace,SFMono-Regular,Menlo,Consolas,monospace);font-weight:700;color:var(--cs-heading,#0f172a);font-size:.9rem;flex:none;letter-spacing:-.01em}
.card-desc{color:var(--cs-muted,#64748b);font-size:.8rem;margin-left:.15rem}
.card-body{padding:.4rem .85rem .9rem}
.seg-adp{background:#99e6da}.seg-capture{background:#a8d5f5}.seg-mask{background:#e3e8ef}.seg-inline{background:#f7c94a}.seg-polytail{background:#f7b1b1}.seg-back{background:#cdc4ee}.seg-insert{background:#f59e0b}
.card-actions{display:flex;gap:.5rem;flex-wrap:wrap;align-items:center;margin:.5rem 0 .2rem}
.copy-btn{border:1px solid var(--cs-border-strong,#cbd5e1);background:var(--cs-bg,#fff);color:var(--cs-text,#1f2937);padding:.34rem .72rem;border-radius:var(--cs-radius-sm,8px);cursor:pointer;font-size:.8rem;font-weight:500;transition:border-color .12s,background .12s,color .12s}
.copy-btn:hover{background:var(--cs-primary-ghost,rgba(13,148,136,.08));border-color:var(--cs-primary,#0d9488);color:var(--cs-primary-strong,#0f766e)}
.copy-btn:focus-visible{outline:2px solid var(--cs-primary,#0d9488);outline-offset:2px}
.card-points{margin:.4rem 0 0;padding-left:1.15rem;color:var(--cs-text,#1f2937);font-size:.85rem}
.card-points li{margin:.16rem 0}
.sg-struct{margin:.45rem 0 .15rem;border:1px solid var(--cs-border,#e2e8f0);border-radius:var(--cs-radius,10px);padding:.55rem .7rem;background:var(--cs-bg-soft,#f8fafc)}
.sg-caption{font-size:.8rem;color:var(--cs-muted,#64748b);margin:.05rem 0 .4rem;font-weight:600;letter-spacing:.01em}
.sg-dl{display:inline-block;border:1px solid var(--cs-primary,#0d9488);background:var(--cs-primary-ghost,rgba(13,148,136,.1));color:var(--cs-primary-strong,#0f766e);padding:.34rem .72rem;border-radius:var(--cs-radius-sm,8px);text-decoration:none;font-size:.8rem;font-weight:600;transition:background .12s,color .12s}
.sg-dl:hover{background:var(--cs-primary-ghost-2,rgba(13,148,136,.18));color:var(--cs-primary-deep,#115e59);text-decoration:none}
.sg-map-svg{width:100%;overflow-x:auto;border:1px solid var(--cs-border,#e2e8f0);border-radius:var(--cs-radius-sm,8px);background:var(--cs-bg,#fff);padding:5px 0}
.sg-map-svg svg{width:100%;min-width:560px;height:auto;display:block}
.feat-band{cursor:pointer;outline:none}
.feat-band rect{transition:stroke-width .12s,filter .12s}
.feat-band:hover rect,.feat-band:focus rect{stroke:var(--cs-primary-deep,#115e59)!important;stroke-width:1.8;filter:brightness(1.05)}
.feat-band.copied rect{stroke:var(--cs-primary,#0d9488)!important;stroke-width:2.4}
.feat-band.feat-struct rect{filter:grayscale(.3)}
.feat-label{font-size:10.5px;fill:#134e4a;font-weight:700;pointer-events:none;paint-order:stroke;stroke:var(--cs-bg,#fff);stroke-width:2.5px;stroke-linejoin:round}
.sg-scale{font-size:11px;font-weight:600;fill:#475569}
.sg-env{font-size:11px;fill:#94a3b8;font-weight:600}
.sg-tip{position:fixed;z-index:80;max-width:340px;background:#0f2e2b;color:#e2e8f0;border:1px solid #155e56;border-radius:var(--cs-radius-sm,8px);padding:9px 11px;font-size:12px;line-height:1.45;pointer-events:none;display:none;font-family:var(--cs-font-sans,system-ui,sans-serif);box-shadow:var(--cs-shadow-lg,0 16px 44px rgba(0,0,0,.4))}
.sg-tip.show{display:block}
.sg-tip h4{margin:0 0 2px;font-size:13px;color:#fff;font-weight:700}
.sg-tip .tip-meta{color:#9fd6cc;font-size:11px;margin:2px 0}
.sg-tip .tip-seq{font-family:var(--cs-font-mono,ui-monospace,Consolas,monospace);font-size:11px;word-break:break-all;background:#134642;border-radius:5px;padding:4px 6px;margin-top:4px;max-height:110px;overflow:auto;color:#e2e8f0}
.card-details{margin:.5rem 0 0;font-size:.82rem}
.card-details>summary{cursor:pointer;color:var(--cs-primary-strong,#0f766e);font-weight:600;font-size:.8rem;list-style:none}
.card-details>summary:hover{color:var(--cs-primary-deep,#115e59)}
.card-details>summary::-webkit-details-marker{display:none}
.sg-table{width:100%;border-collapse:collapse;font-size:.78rem;margin-top:.6rem}
.sg-table th,.sg-table td{border-bottom:1px solid var(--cs-border,#e2e8f0);padding:.3rem .45rem;text-align:left;vertical-align:top}
.sg-table th{color:var(--cs-muted,#64748b);font-weight:600;background:var(--cs-table-head,#f1f5f9)}
.sg-table tbody tr:hover{background:var(--cs-primary-ghost,rgba(13,148,136,.04))}
.sg-role{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:5px;vertical-align:middle}
"""

_JS = """
(function(){
  var D = JSON.parse(document.getElementById('adapdata').textContent);
  var host = document.getElementById('adapp');
  var box = document.getElementById('adsearch');
  var countEl = document.getElementById('adcount');

  function esc(s){ return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  var ROLE_COLOR = {adp:'#99e6da', capture:'#a8d5f5', mask:'#e3e8ef', inline:'#f7c94a', polytail:'#f7b1b1', back:'#cdc4ee', insert:'#f59e0b'};
  var ROLE_LBL = {adp:'adapter', capture:'UMI', mask:'mask', inline:'barcode', polytail:'poly-tail', back:'read-through', insert:'insert'};
  // Short, readable band label from a feature name (drop redundant prefixes /
  // flowcell/index scaffolding) so the diagram isn't a wall of generic 'ADAPTER'.
  function shortName(f){
    var n = f.name || ROLE_LBL[f.role] || f.role;
    n = String(n).replace(/library |TruSeq sRNA |TruSeq |Nextera |5&prime;|<i>/gi, '').replace(/ \(read\)/gi, '').replace(/\s*\(\d+N?\).*$/i, '').trim();
    if (f.role === 'insert'){ return 'INSERT'; }
    if (f.role === 'mask'){ return 'MASK'; }
    if (f.role === 'capture'){ return 'UMI'; }
    if (f.role === 'inline'){ return 'BC'; }
    if (f.role === 'polytail'){ return 'poly'; }
    return n.length ? n.toUpperCase() : ROLE_LBL[f.role].toUpperCase();
  }
  // Flowcell anchors / index placeholders are library-prep, NOT read-visible.
  function isStructural(f){
    var n = String(f.name || '');
    return /flowcell/i.test(n) || /index/i.test(n) || /^P\d/i.test(n.trim());
  }

  // One shared, floating detail tooltip for the whole viewer. Delegated over the
  // SVG feature bands so a single element serves every card.
  var tip = document.createElement('div');
  tip.className = 'sg-tip';
  document.body.appendChild(tip);
  function featHtml(f){
    return 'data-role="'+esc(f.role)+'" data-name="'+esc(f.name)+'" data-start="'+f.start+'" data-end="'+f.end+'" data-len="'+f.len+'" data-side="'+esc(f.side||'')+'" data-seq="'+esc(f.role==='insert' ? ('N'+f.len) : f.seq)+'"';
  }
  function tipText(f){
    var loc = f.start + '..' + f.end;
    var role = ROLE_LBL[f.role] || f.role;
    return '<h4>'+esc(f.name)+'</h4>'
      + '<div class="tip-meta">'+esc(role)+' &middot; '+loc+' &middot; '+f.len+' bp'+(f.side ? ' &middot; '+esc(f.side) : '')+'</div>'
      + '<div class="tip-seq">'+(f.role==='insert' ? 'N'+f.len : esc(f.seq))+'</div>';
  }
  // Build a band-feature surrogate from the element's data-* attributes so the
  // delegated handlers don't need the original object handy.
  function bandFeat(el){
    return {
      role: el.getAttribute('data-role'),
      name: el.getAttribute('data-name'),
      start: el.getAttribute('data-start'),
      end: el.getAttribute('data-end'),
      len: el.getAttribute('data-len'),
      side: el.getAttribute('data-side'),
      seq: el.getAttribute('data-seq')
    };
  }
  function showTip(el, evt){
    tip.innerHTML = tipText(bandFeat(el));
    tip.classList.add('show');
    // keep the tooltip reasonably within the viewport
    var r = el.getBoundingClientRect();
    var tTop = r.bottom + 4;
    var h = tip.offsetHeight;
    if (tTop + h > window.innerHeight - 8){ tTop = r.top - h - 4; }
    tip.style.left = Math.max(8, Math.min(r.left, window.innerWidth - tip.offsetWidth - 8)) + 'px';
    tip.style.top = tTop + 'px';
  }
  function hideTip(){ tip.classList.remove('show'); }

  function mapSvg(c){
    var k = c.construct, size = k.size, feats = k.features;
    var w = 920, h = 70, x0 = 8, x1 = w - 8;
    function sx(p){ return x0 + (x1 - x0) * (p - 1) / Math.max(1, size - 1); }
    var s = '<div class="sg-map-svg"><svg viewBox="0 0 '+w+' '+h+'" width="100%" role="list" aria-label="Structure diagram">';
    s += '<line x1="'+x0+'" y1="34" x2="'+x1+'" y2="34" stroke="#cbd5e1" stroke-width="16" stroke-linecap="round"/>';
    s += '<text class="sg-env" x="'+x0+'" y="64">5&prime;</text>';
    s += '<text class="sg-env" x="'+x1+'" y="64" text-anchor="end">3&prime;</text>';
    for (var i=0;i<feats.length;i++){
      var f = feats[i];
      var a = sx(f.start), b = sx(f.end);
      if (b - a < 6){ var mid=(a+b)/2; a=mid-3; b=mid+3; }
      var col = ROLE_COLOR[f.role] || '#e0e0e0';
      // Flowcell/index scaffolding is library-prep, NOT read-visible: dim it so the
      // read-visible adapter / UMI / mask / insert structure is the focus.
      var struct = isStructural(f);
      var bandCls = struct ? 'feat-band feat-struct' : 'feat-band';
      var op = struct ? ' opacity="0.36"' : '';
      // Each band is one feature: hovering/focusing shows its detail, clicking copies its sequence.
      s += '<g class="'+bandCls+'" tabindex="0" role="listitem" data-copy="'+esc(f.role==='insert' ? ('N'+f.len) : f.seq)+'" '+featHtml(f)+'>';
      s += '<rect x="'+a+'" y="27" width="'+(b-a)+'" height="14" rx="3" fill="'+col+'" stroke="#94a3b8" stroke-width="1"'+op+'/>';
      if (f.strand === '-'){ s += '<path d="M '+(a+8)+' 27 L '+(a)+' 34 L '+(a+8)+' 41 z" fill="'+col+'" stroke="#64748b"'+op+'/>'; }
      else { s += '<path d="M '+(b-8)+' 27 L '+(b)+' 34 L '+(b-8)+' 41 z" fill="'+col+'" stroke="#64748b"'+op+'/>'; }
      // Inline specific label when the read-visible band can carry it — only when
      // it actually fits, so short names never get crushed or overlap.
      if (!struct){ var lbl = shortName(f); var labW = lbl.length*6.6 + 12; if ((b - a) >= labW){ s += '<text class="feat-label" x="'+((a+b)/2)+'" y="38" text-anchor="middle">'+esc(lbl)+'</text>'; } }
      s += '</g>';
    }
    s += '<text class="sg-scale" x="'+(x0+x1)/2+'" y="63" text-anchor="middle">'+size+' bp &middot; linear</text>';
    s += '</svg></div>';
    return s;
  }

  function caption(c){
    var k = c.construct;
    var ins = 0;
    for (var i=0;i<k.features.length;i++){ if (k.features[i].role === 'insert'){ ins++; } }
    return '<div class="sg-caption">'+k.size+' bp &middot; 5&prime;&rarr;3&prime; top strand' + (ins ? ' &middot; '+ins+' insert placeholder' : '') + '</div>';
  }

  function featTable(c){
    var feats = c.construct.features;
    var s = '<table class="sg-table"><thead><tr><th>Feature</th><th>Location</th><th>Size</th><th>Role</th><th>Arm</th><th>Sequence</th></tr></thead><tbody>';
    if (!feats.length){ s += '<tr><td colspan=6>no features</td></tr>'; }
    for (var i=0;i<feats.length;i++){
      var f = feats[i];
      var sq = f.seq;
      if (f.role === 'insert'){ sq = 'N\u00d7'+f.len; }
      s += '<tr><td>' + esc(f.name) + '</td><td>' + f.start + '..' + f.end + '</td><td>' + f.len + '</td><td><span class="sg-role" style="background:'+(ROLE_COLOR[f.role]||'#e0e0e0')+'"></span>' + ROLE_LBL[f.role] + '</td><td>' + (f.side||'&mdash;') + '</td><td style="font-family:monospace;max-width:220px;word-break:break-all">' + esc(sq.length > 40 ? sq.slice(0,40)+'…' : sq) + '</td></tr>';
    }
    s += '</tbody></table>';
    return s;
  }

  function structBlock(c){
    var k = c.construct;
    if (!k || !k.size){ return ''; }
    // One primary interactive diagram; the per-base sequence is gone and the full
    // feature table folds behind a small collapse, so a card reads as a single
    // clear construct rather than three redundant stacked representations.
    return '<div class="sg-struct">' + caption(c) + mapSvg(c)
      + '<details class="card-details"><summary>Feature details</summary>' + featTable(c) + '</details>'
      + '</div>';
  }

  function card(c){
    var alias = c.alias ? ' data-alias="1"' : '';
    var cmd = 'cutseq -A ' + c.name;
    var points = '';
    if (c.points && c.points.length){ points = '<ul class="card-points">' + c.points.map(function(p){return '<li>'+esc(p)+'</li>';}).join('') + '</ul>'; }
    var s = '';
    s += '<details class="adapter-card"'+alias+'>';
    s += '<summary data-copy="'+esc(c.scheme)+'"><span class="card-name">'+esc(c.name)+'</span><span class="card-desc">'+esc(c.desc)+'</span>'+(c.alias?'<span class="card-desc">[alias of '+esc(c.alias)+']</span>':'')+'</summary>';
    s += '<div class="card-body">';
    // One primary interactive structure diagram (mapSvg) + the hidden-by-default
    // feature table inside structBlock. The old colored seg-row chip list is gone.
    s += structBlock(c);
    s += '<div class="card-actions">';
    s += '<button class="copy-btn" data-copy="'+esc(c.scheme)+'">copy scheme</button>';
    s += '<button class="copy-btn" data-copy="'+esc(cmd)+'">copy: '+esc(cmd)+'</button>';
    s += '<a class="sg-dl" href="constructs/'+encodeURIComponent(c.name)+'.gb" download>Download .gb (GenBank)</a>';
    s += '</div>' + points + '</div></details>';
    return s;
  }

  function render(){
    var q = (box.value || '').trim().toLowerCase();
    var n = 0;
    host.innerHTML = D.schemes.filter(function(c){
      if (!q){ return true; }
      var hay = (c.name + ' ' + c.desc + ' ' + c.scheme).toLowerCase();
      if (hay.indexOf(q) !== -1){ return true; }
      for (var i=0;i<c.tokens.length;i++){ if (c.tokens[i].s && c.tokens[i].s.toLowerCase().indexOf(q) !== -1){ return true; } }
      return false;
    }).map(function(c){ n++; return card(c); }).join('');
    countEl.textContent = n + ' / ' + D.schemes.length + ' schemes';
  }
  box.addEventListener('input', render);
  render();

  // Tooltip on hover / focus over a feature band (delegated).
  function track(e){
    var el = e.target.closest ? e.target.closest('.feat-band') : null;
    if (el){ showTip(el, e); } else { hideTip(); }
  }
  document.addEventListener('mousemove', track);
  document.addEventListener('mouseleave', hideTip);
  document.addEventListener('focusin', function(e){
    var el = e.target.closest ? e.target.closest('.feat-band') : null;
    if (el){ showTip(el, e); }
  }, true);
  document.addEventListener('focusout', hideTip, true);

  document.addEventListener('click', function(e){
    var el = e.target.closest ? e.target.closest('[data-copy]') : null;
    if (!el){ return; }
    var val = el.getAttribute('data-copy');
    function done(){
      // Only the 'copy scheme'/'copy: cmd' BUTTONS are action labels whose
      // text can be briefly replaced by 'copied!'. A colored segment's text is
      // the actual scheme sequence — never overwrite it; flash it instead.
      var isBtn = el.classList && el.classList.contains('copy-btn');
      if (isBtn){
        var old = el.getAttribute('data-label');
        if (!old){ el.setAttribute('data-label', el.textContent); el.textContent = 'copied!'; setTimeout(function(){ el.textContent = el.getAttribute('data-label'); el.removeAttribute('data-label'); }, 900); }
      } else {
        var old = el.getAttribute('data-fb');
        if (!old){ el.setAttribute('data-fb','1'); el.classList.add('copied'); setTimeout(function(){ el.classList.remove('copied'); el.removeAttribute('data-fb'); }, 700); }
      }
    }
    if (navigator.clipboard && navigator.clipboard.writeText){ navigator.clipboard.writeText(val).then(done, function(){ legacyCopy(); }); }
    else { legacyCopy(); }
    function legacyCopy(){
      var t = document.createElement('textarea'); t.value = val; t.setAttribute('readonly','');
      t.style.position = 'fixed'; t.style.left = '-9999px';
      document.body.appendChild(t); t.select();
      try{ document.execCommand('copy'); }catch(err){}
      document.body.removeChild(t); done();
    }
    e.stopPropagation();
  });
})();
"""


def _genbank(name, desc, scheme, seq, feats):
    """Render a GenBank (.gb) record for a schematic library construct."""
    size = len(seq)
    L = []
    L.append(f"LOCUS       {name:<14}{size:>10} bp    DNA     linear   {_gb_date():>11} SYN")
    L.append(f"// cutseq built-in scheme: {name} ({desc})")
    L.append(f"// grammar: {scheme}")
    L.append("FEATURES             Location/Qualifiers")
    L.append(f"     source          1..{size}")
    L.append('                     /organism="synthetic construct"')
    L.append('                     /mol_type="genomic DNA"')
    L.append('                     /note="cutseq '+ (name or "") + ' library scheme; adapter/barcode/UMI structure"')
    for f in feats:
        loc = f"{f['start']}..{f['end']}"
        L.append(f"     {f['type']:<16}{loc}")
        L.append(f'                     /label="{_gb_quote(f["name"])}"')
        note = f"cutseq role: {f['role']}"
        if f.get("side"):
            note += f"; arm: {f['side']}"
        if f.get("strand") == "-":
            note += "; 3' read-through"
        L.append(f'                     /note="{_gb_quote(note)}"')
    L.append("ORIGIN")
    s = seq.upper()
    for i in range(0, size, 60):
        chunk = s[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        L.append(f"{i:>9} {groups}")
    L.append("//")
    return "\n".join(L) + "\n"


def _gb_quote(s):
    return s.replace('"', "'").replace(chr(10), " ")


def _gb_date():
    """A deterministic GenBank LOCUS date.

    Using a fixed date keeps the generated ``.gb`` byte-identical across runs and
    machines (local vs CI), so the ``update-adapter-docs`` pre-commit hook never
    reports a spurious change. The LOCUS date is display-only metadata.
    """
    return "01-JAN-1970"


def _write_exports(name, construct, constructs_dir):
    """Write a GenBank (.gb) file for one scheme."""
    constructs_dir.mkdir(parents=True, exist_ok=True)
    gb = _genbank(name, "", "", construct["seq"], construct["features"])
    (constructs_dir / f"{name}.gb").write_text(gb, encoding="utf-8")



_DARK_CSS = """
@media (prefers-color-scheme: dark) {
  /* The viewer mostly inherits --cs-* tokens from the shared shell, so only the
     literal SVG fills / tooltip / tiny accents need explicit handling here. */
  .ad-toolbar{background:#101824;border-color:#26303f}
  .ad-search{background:#0b1017;color:#e2e8f0;border-color:#3b475e}
  .ad-search::placeholder{color:#64748b}
  .ad-count{color:#94a3b8}
  .ad-legend{background:#0e1526;border-color:#26303f}
  .adleg{color:#cbd5e1}
  .adapter-card{background:#0b1017;border-color:#26303f}
  .adapter-card>summary{background:#0e1526}
  .adapter-card>summary:hover{background:#12333a}
  .card-name{color:#e2e8f0}
  .card-desc{color:#94a3b8}
  .card-points{color:#cbd5e1}
  .sg-struct{background:#0e1526;border-color:#26303f}
  .sg-caption{color:#9fb2c3}
  .sg-map-svg{background:#0b1017;border-color:#26303f}
  .sg-table th{background:#131b2a;color:#cbd5e1}
  .sg-table th,.sg-table td{border-color:#26303f}
  .sg-table tbody tr:hover{background:#12333a}
  .card-details>summary{color:#5eead4}
  .copy-btn{background:#131b2a;color:#e2e8f0;border-color:#3b475e}
  .copy-btn:hover{background:#12333a;border-color:#5eead4;color:#e2e8f0}
  .sg-dl{background:rgba(45,212,191,.12);color:#5eead4;border-color:#14b8a6}
  .sg-dl:hover{background:rgba(45,212,191,.2);color:#99f6e4}
  .feat-label{fill:#eafcf9;stroke:#0b1017}
  .sg-scale{fill:#9fb2c3}
  .sg-env{fill:#64748b}
  .sg-tip{background:#0f2e2b}
}
"""

# --- standalone viewer page (self-contained, Jekyll/Liquid-safe) -------------
_INTRO = (
    "<p>CutSeq bundles built-in library schemes &mdash; list them in a terminal with "
    "<code>cutseq --list-adapters</code>, and use a name (or a custom grammar string) "
    "via <code>-A/--adapter-scheme</code> (<code>-a</code> is an alias). This page is "
    "a live viewer: <strong>search</strong> to filter, <strong>hover</strong> a colored "
    "segment to see its role, <strong>click</strong> any segment or button to copy it.</p>"
)

_NOTES = (
    "<h2>Inline-barcode auto-detection</h2>\n"
    "<p>Inline barcodes are written in lowercase (<code>acgt...</code>). If you write "
    "one in uppercase, it merges with the adjacent sequencing primer into one adapter "
    "run. CutSeq checks the two outermost adapters of a scheme against a curated "
    "database of Illumina / BGI (MGI) sequencing primers "
    "(<code>cutseq --list-primers</code>) and reclassifies any fixed uppercase "
    "sequence adjacent to (or between) recognized primers as an inline barcode, with "
    "a warning. Custom schemes without known primers are never altered. Disable with "
    "<code>--no-auto-inline</code>.</p>\n"
    "<h2>Discarding reads with a reason tag</h2>\n"
    "<p>A single <code>-d</code> output collects every rejected read with the reason "
    "in the name (<code>reason=too_short</code>, <code>reason=too_many_n</code>, "
    "<code>reason=low_quality</code>, <code>reason=no_barcode</code>), gated by "
    "<code>-m/--min-length</code>, <code>--max-n</code>, <code>--min-avg-quality</code> "
    "and <code>--ensure-inline-barcode</code>. When no <code>-d</code>/<code>-O</code> "
    "is given, discard files are auto-named from the inputs.</p>"
)

# Folded from docs/spatial.html (deleted): the DBiT-seq barcode-arm library walkthrough.
# Kept inside _NOTES so a generator re-run reproduces adapters.html byte-for-byte.
# Folded from docs/spatial.html (deleted): the DBiT-seq barcode-arm library walkthrough.
# Kept inside _NOTES so a generator re-run reproduces adapters.html byte-for-byte.
_DBIT_SECTION = r"""
<section id="spatial-dbit-seq">
<h2 id="parsing-spatial-libraries">Parsing spatial (DBiT-seq) libraries</h2>
<p>This page is the introduction to handling <strong>spatial sequencing</strong> data with CutSeq: what the library construct looks like, which read carries what, the one-line command that trims everything and captures the spatial barcodes into the read names, and how to interpret the QC that results from it.</p>
<p>It is written for the <strong>DBiT-seq</strong> scheme (<code>cutseq -A DBITSEQ</code>, spatial RNA, 50&times;50 barcode grid), but the same walk-through applies to any barcode-arm library expressed as a CutSeq scheme.</p>
<hr>
<h3 id="1-the-library-construct">1. The library construct</h3>
<p>DBiT-seq puts a row of fixed-length oligos — the <strong>barcode arm</strong> — on one side of each molecule. The full construct, top strand 5&prime;&rarr;3&prime;:</p>
<pre><code class="language-plaintext">P5(29) | i5(8) | R1-primer(34) | TSO(23) | insert | polyA/T | UMI(12)
       | rc(L1)(30) | BarcodeA(8) | rc(L2)(30) | BarcodeB(8)
       | rc(handle)(22) | rc(R2-primer)(34) | i7(8) | rc(P7)(24)</code></pre>
<p>The sequencing primers just upstream of each read are <strong>never part of the reads</strong> — sequences start downstream of the priming site. So:</p>
<table><thead><tr><th>Read</th><th>What it actually sequences</th></tr></thead><tbody><tr><td><strong>R1</strong> (from the P5 side)</td><td><strong>TSO</strong> (23 nt template-switch oligo) + cDNA <strong>insert</strong> (variable length), often with a 5&prime; G-stretch right after the TSO</td></tr><tr><td><strong>R2</strong> (from the P7 side)</td><td><strong>barcode arm</strong> then insert read-through: handle(22) → <strong>BarcodeB</strong>(8) → linker2(30) → <strong>BarcodeA</strong>(8) → linker1(30) → <strong>UMI</strong>(12) → polyA/T → 3&prime; end of the insert</td></tr></tbody></table>
<div class="callout note"><div class="callout-title">Note</div><p>The insert is on R1 (5&prime;&rarr;3&prime;, top strand). R2 only sees the <strong>3&prime; tail</strong> of the insert that reads through past the arm — and it is capped by the sequencer read length (see <a href="#5-interpreting-the-read-length-distribution">§5 – interpreting read lengths</a>).</p></div>
<h3 id="2-one-command-to-parse-everything">2. One command to parse everything</h3>
<pre><code class="language-bash">cutseq -A DBITSEQ \
  -R '{id}_BCB:{3}_BCA:{2}_UMI:{1}' \
  -O mysample \
  mysample_R1.fastq.gz mysample_R2.fastq.gz</code></pre>
<p>Fast, single-threaded by default (<code>-t N</code> to parallelize). That one line:</p>
<ul><li>trims the <strong>TSO + 5&prime; G-stretch off R1</strong>,</li><li>walks the <strong>barcode arm on R2</strong> (handle, both linkers, poly tail),</li><li><strong>captures BarcodeB, BarcodeA and the UMI</strong> and writes them into every read name, e.g.:</li></ul>
<pre><code class="language-plaintext">@LH00699:99:2533WTLT4:6:1101:36114:1055_BCB:AGAGTCAA_BCA:ACCTCCAA_UMI:TGCACACCAACC</code></pre>
<h4 id="what-the-built-in-dbitseq-scheme-means">What the built-in <code>DBITSEQ</code> scheme means</h4>
<pre><code class="language-plaintext">AAGCAGTGGTATCAACGCAGAGTGAATGGG...GGG : N12 AGTCGTACGCCGATGCGAAACATCGGCCAC
                                            N8  CGAATGCTCTGGCCTCTCAAGCACGTGGAT
                                            N8  AGATGCGAGAAGCCAACGCTTG</code></pre>
<table><thead><tr><th>Token</th><th>Side</th><th>Meaning in CutSeq</th></tr></thead><tbody><tr><td><code>AAGCAGTGGTATCAACGCAGAGT</code></td><td>left (<code>:</code>⇒R1, as-is)</td><td>TSO — R1&rsquo;s real read-5&prime; scaffold, trimmed off</td></tr><tr><td><code>GAATGGG...GGG</code></td><td>R1</td><td>auto-detected <strong>5&prime; poly-Run</strong> (G-stretch / GAATGGG) leftmost-anchored — trims the template-switch remnant; <code>G&lt;k&gt;</code> tunes the minimum run</td></tr><tr><td><code>N12</code></td><td>right (matched on R2)</td><td><strong>UMI</strong> capture &rarr; <code>{1}</code></td></tr><tr><td>first 30-nt run</td><td>right</td><td>linker 1, trimmed</td></tr><tr><td><code>N8</code></td><td>right</td><td><strong>BarcodeA</strong> capture &rarr; <code>{2}</code></td></tr><tr><td>second 30-nt run</td><td>right</td><td>linker 2, trimmed</td></tr><tr><td><code>N8</code></td><td>right</td><td><strong>BarcodeB</strong> capture &rarr; <code>{3}</code></td></tr><tr><td>final 22-nt run</td><td>right</td><td>handle (rc form), trimmed</td></tr></tbody></table>
<p>The <code>-R '{id}_BCB:{3}_BCA:{2}_UMI:{1}'</code> template is read-name renaming — see the <a href="rename.html">rename page</a> for the full syntax (transforms, <code>rc()</code>, per-side captures, etc.).</p>
<h3 id="3-outputs">3. Outputs</h3>
<p>With <code>-O mysample</code> you get (input file names are used if <code>-O</code> is omitted):</p>
<pre><code class="language-plaintext">mysample_trimmed_R1.fastq.gz   # insert read, barcodes in header, TSO/G-stretch removed
mysample_trimmed_R2.fastq.gz   # barcode arm trimmed, barcodes in header, UMI in header
mysample_discard_R1.fastq.gz   # rejected reads (both mates kept together)
mysample_discard_R2.fastq.gz</code></pre>
<p>Rejected reads carry a <code>reason=</code> tag in the name: <code>too_short</code> (below <code>-m/--min-length</code>, default 20), <code>too_many_n</code> (<code>--max-n</code>), <code>low_quality</code> (<code>--min-avg-quality</code>), <code>no_barcode</code> (only with <code>--ensure-inline-barcode</code>). Add <code>--json-file mysample.json</code> for a full trimming report with counts/statistics per step.</p>
<div class="callout note"><div class="callout-title">Note</div><p>Because trimming is performed on paired, in-order reads, R1 (insert) and R2 (barcode) mates stay in sync — the barcode in an R2 header belongs to the insert sequence in the same R1 record.</p></div>
<h3 id="4-going-from-trimmed-reads-to-a-spatial-expression-matrix">4. Going from trimmed reads to a spatial expression matrix</h3>
<ol><li><strong>Parse the header</strong> — split each read name on <code>_BCB:</code>, <code>_BCA:</code>, <code>_UMI:</code> and read the 8+8+12 nt.</li><li><strong>Map (BarcodeB, BarcodeA) → pixel</strong> with the design&rsquo;s barcode whitelist (the two 8-nt sets used to build the 50&times;50 grid; allow 1–2 nt Levenshtein distance if the protocol tolerates it).</li><li><strong>Align R1</strong> to the reference and assign transcripts/features per pixel (m6A-ARTR-style methods typically report per-feature, per-pixel read counts).</li><li><strong>Deduplicate by UMI</strong> — collapse reads sharing the same pixel + feature (+ alignment position) and the same 12-nt UMI.</li></ol>
<p>A minimal header parser / demultiplexer (paired trimmed R1+R2, barcodes in the R1 name):</p>
<pre><code class="language-python">import re, gzip

PAT = re.compile(r"_BCB:([ACGT]{8})_BCA:([ACGT]{8})_UMI:([ACGT]{12})")

def records(fq):
    with gzip.open(fq, "rt") as f:
        lines = [ln.strip() for ln in f]
    for i in range(0, len(lines), 4):
        yield lines[i], lines[i + 1]                # (name, sequence)

with open("pixels.tsv", "w") as out:
    for (r1n, r1s), (r2n, r2s) in zip(
            records("mysample_trimmed_R1.fastq.gz"),
            records("mysample_trimmed_R2.fastq.gz")):
        m = PAT.search(r1n)
        if not m:
            continue
        bcb, bca, umi = m.groups()
        out.write(f"{r1n.split()[0]}\t{bcb}\t{bca}\t{umi}\t{r1s}\n")</code></pre>
<h3 id="5-interpreting-the-read-length-distribution">5. Interpreting the read-length distribution (the R2 “41-nt peak”)</h3>
<p>After trimming, R1 lengths reflect the real insert-length spread (here ~9–123 nt, peaking at ~11–15). <strong>R2 is different</strong>: the fixed barcode arm occupies 110 nt (handle 22 + BCB 8 + linker2 30 + BCA 8 + linker1 30 + UMI 12) of the sequenced read, so the arm + UMI + read-through has to fit in the read.</p>
<p>For a 151-cycle run:</p>
<pre><code class="language-plaintext">151 (read length) − 110 (arm: handle22 + BCB8 + L2 30 + BCA8 + L1 30 + UMI12) = 41 nt</code></pre>
<p>Every molecule whose insert is <strong>longer than ~41 nt on the R2 side</strong> therefore hits the physical end of the read and piles up at exactly <strong>41 nt</strong> after trimming (<code>raw[-41:]</code> in 99% of cases). In practice this can be ~18% of all R2 reads — it is a <strong>read-truncation artefact, not a biological 41-nt fragment</strong>: the sequences are fully diverse and their R1 mates are long.</p>
<div class="callout warning"><div class="callout-title">Warning</div><p>Do <strong>not</strong> read R2&rsquo;s length axis as fragment length. R2 is capped by the sequencer, not by the insert. If you need the full 3&prime; R2-side sequence read length 151→250 (<code>--r2-primer</code> is informational only), or rely on <strong>R1</strong> for insert biology.</p></div>
<p>Correspondingly, in a per-position <strong>base-composition &times; depth</strong> pileup:</p>
<ul><li><strong>R1</strong> — a clean, slightly composition-skewed stack that runs out around the longest insert lengths;</li><li><strong>R2</strong> — strong constant composition across the arm-bearing positions, then a tail of read-through cDNA; expect the depth to collapse once reads that reached the 151-nt end stop contributing.</li></ul>
<p>The built-in defaults are gentle (<code>-m 20</code>, <code>-q 20</code>, no <code>--max-n</code> / <code>--min-avg-quality</code> by default); set them from your QC thresholds and check the <code>reason=</code> distribution in the discard files.</p>
<h3 id="6-notes-extras">6. Notes &amp; extras</h3>
<ul><li>Legacy name: <code>-A M6AARTR</code> still works as an alias for <code>-A DBITSEQ</code> — prefer <code>-A DBITSEQ</code>.</li><li>The same library as a fully inline custom scheme (no built-in needed): <code>cutseq -A "AAGCAGTGGTATCAACGCAGAGTGAATGGG...GGG:N12AGTCGTACGCCGATGCGAAACATCGGCCACN8CGAATGCTCTGGCCTCTCAAGCACGTGGATN8AGATGCGAGAAGCCAACGCTTG" -R '{id}_BCB:{3}_BCA:{2}_UMI:{1}' ...</code></li><li>Dry-run inspection: <code>cutseq -A DBITSEQ --graph-vertical R1.fq.gz R2.fq.gz</code> prints the parsed trimming steps without writing any output — useful to confirm the arm layout before running the real job.</li><li><code>--r1-primer</code> / <code>--r2-primer</code> document the sequencing primers but are never trimmed from reads (reads start downstream of them).</li><li>See <a href="quickstart.html">Quick Start</a>, <a href="adapters.html">Adapter schemes</a> and <a href="rename.html">Read-name renaming</a> for more.</li></ul>
</section>
"""

_NOTES = _NOTES + _DBIT_SECTION
_PAGE_TEMPLATE = (
    "<!DOCTYPE html>\n"
    "<html lang=\"en\">\n"
    "<head>\n"
    "<meta charset=\"utf-8\">\n"
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
    "<title>Adapter Schemes &mdash; CutSeq</title>\n"
    "<meta name=\"description\" content=\"All built-in CutSeq library schemes, with an interactive viewer.\">\n"
    "<link rel=\"stylesheet\" href=\"assets/css/site.css\">\n"
    "<link rel=\"icon\" href=\"favicon.ico\">\n"
    "<script async src=\"https://www.googletagmanager.com/gtag/js?id=G-EFHY8WFMHK\"></script>\n"
    "<script>\n"
    "  window.dataLayer = window.dataLayer || [];\n"
    "  function gtag(){{dataLayer.push(arguments);}}\n"
    "  gtag('js', new Date());\n"
    "  gtag('config', 'G-EFHY8WFMHK', {{ 'anonymize_ip': true }});\n"
    "</script>\n"
    "<style>{css}</style>\n"
    "<style>{dark_css}</style>\n"
    "</head>\n"
    "<body>\n"
    "<header class=\"site-header\">\n"
    "  <div class=\"site-title-box\"><a class=\"site-title\" href=\"index.html\"><img class=\"site-logo\" src=\"logo.svg\" alt=\"cutseq logo\"><span>cutseq</span></a></div>\n"
    "  <div class=\"site-aux\"><a href=\"https://github.com/y9c/cutseq\" target=\"_blank\" rel=\"noopener\">View on GitHub</a></div>\n"
    "</header>\n"
    "<div class=\"site-shell\">\n"
    "  <nav class=\"site-nav\" aria-label=\"Main\">\n"
    "    <ul class=\"nav-list\">\n"
    "      <li class=\"nav-item\"><a href=\"index.html\">Home</a></li>\n"
    "      <li class=\"nav-item\"><a href=\"quickstart.html\">Quick Start</a></li>\n"
    "      <li class=\"nav-item active\"><a aria-current=\"page\" href=\"adapters.html\">Adapter Schemes</a></li>\n"
    "      <li class=\"nav-item\"><a href=\"rename.html\">Read Name Renaming</a></li>\n"
    "    </ul>\n"
    "  </nav>\n"
    "  <div class=\"site-main\">\n"
    "    <main class=\"site-content site-content-wide\" id=\"content\">\n"
    "<h1>Adapter Schemes</h1>\n"
    "{intro}\n"
    "{legend}\n"
    "<div class=\"ad-toolbar\">\n"
    "<input id=\"adsearch\" class=\"ad-search\" type=\"search\" "
    "placeholder=\"Search schemes\u2026\" />\n"
    "<span id=\"adcount\" class=\"ad-count\">{count} schemes</span>\n"
    "</div>\n"
    "<div id=\"adapp\"></div>\n"
    "<script type=\"application/json\" id=\"adapdata\">{data_json}</script>\n"
    "<script>{js}</script>\n"
    "{notes}\n"
    "    </main>\n"
    "  </div><!-- /site-main -->\n"
    "</div><!-- /site-shell -->\n"
    "<footer class=\"site-footer\">\n"
    "  <p>Copyright &copy; 2021-2026, Chang Ye</p>\n"
    "  <a href=\"#content\" class=\"back-to-top\" id=\"back-to-top\">Back to top</a>\n"
    "</footer>\n"
    "<script src=\"assets/js/site.js\"></script>\n"
    "</body>\n"
    "</html>\n"
)


def _page_html(schemes, data_json, css, dark_css, js, legend):
    """Return the complete shared-shell HTML document for the adapters viewer."""
    return _PAGE_TEMPLATE.format(
        css=css,
        dark_css=dark_css,
        intro=_INTRO,
        legend=legend,
        count=len(schemes),
        data_json=data_json,
        js=js,
        notes=_NOTES,
    )


def main():
    schemes = _load_schemes()
    if not schemes:
        sys.exit("Error: no schemes loaded from adapters.toml")

    legend = _legend_html()
    payload = {"schemes": schemes, "roles": {k: v[0] for k, v in ROLES.items()},
               "legend": legend}
    data_json = json.dumps(payload, ensure_ascii=True)

    # write downloadable construct files (.gb) into docs/constructs
    constructs_dir = project_root / "docs" / "constructs"
    for s in schemes:
        _write_exports(s["name"], s["construct"], constructs_dir)
    print(f"wrote constructs for {len(schemes)} schemes to {constructs_dir}")

    out_html = project_root / "docs" / "adapters.html"
    page = _page_html(schemes, data_json, _CSS, _DARK_CSS, _JS, legend)
    out_html.write_text(page, encoding="utf-8")
    print(f"wrote {out_html}: {len(schemes)} schemes")


if __name__ == "__main__":
    main()