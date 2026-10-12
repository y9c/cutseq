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
.ad-toolbar{display:flex;gap:1rem;align-items:center;flex-wrap:wrap;margin:1rem 0;padding:12px 14px;background:#f6f8fa;border:1px solid #e2e8f0;border-radius:10px}
.ad-search{flex:1 1 260px;font-size:15px;padding:8px 12px;border:1px solid #cbd5e1;border-radius:8px;min-width:200px}
.ad-count{font-size:.9rem;color:#475569;white-space:nowrap}
.ad-legend{display:flex;gap:1.1rem;flex-wrap:wrap;align-items:center;margin:.2rem 0 1rem;padding:.6rem 1rem;background:#fff;border:1px solid #e2e8f0;border-radius:10px}
.adleg{display:flex;align-items:center;gap:.4rem;font-size:.82rem;color:#334155}
.leg-swatch{width:13px;height:13px;border-radius:3px;display:inline-block}
.adapter-card{margin:.6rem 0;border:1px solid #e2e8f0;border-radius:10px;background:#fff;overflow:hidden}
.adapter-card>summary{display:flex;gap:.8rem;align-items:baseline;cursor:pointer;padding:.5rem .75rem;list-style:none;flex-wrap:wrap}
.adapter-card>summary::-webkit-details-marker{display:none}
.adapter-card[data-alias]{border-left:4px solid #f59e0b}
.card-name{font-family:monospace;font-weight:700;color:#0f172a;font-size:.95rem;flex:none}
.card-desc{color:#64748b;font-size:.83rem;margin-left:.2rem}
.card-body{padding:.35rem .75rem .75rem}
.seg-row{display:flex;flex-wrap:wrap;align-items:center;gap:3px;font-family:monospace;font-size:13px;padding:.4rem 0}
.seg{cursor:pointer;padding:4px 7px;border-radius:4px;border:1px solid transparent;white-space:nowrap;position:relative}
.seg:hover{border-color:#94a3b8;box-shadow:0 0 0 2px rgba(59,130,246,.25)}
.seg-adp{background:#A8E6CF}.seg-capture{background:#B2EBF2}.seg-mask{background:#DCDCDC}.seg-inline{background:#FFD700}.seg-polytail{background:#FFB3BA}.seg-back{background:#C3B1E1}.seg-insert{background:#FF6F61;color:#fff;font-weight:700}
.card-actions{display:flex;gap:.5rem;flex-wrap:wrap;margin:.3rem 0 .5rem}
.copy-btn{border:1px solid #cbd5e1;background:#f8fafc;color:#0f172a;padding:.32rem .7rem;border-radius:6px;cursor:pointer;font-size:.8rem}
.copy-btn:hover{background:#eef2ff;border-color:#2563eb}
.card-points{margin:.2rem 0 0;padding-left:1.1rem;color:#334155;font-size:.86rem}
.card-points li{margin:.14rem 0}
.tooltiptext{position:absolute;bottom:110%;left:0;background:#0f172a;color:#e2e8f0;padding:3px 7px;border-radius:5px;font-size:11px;white-space:nowrap;display:none;z-index:20;font-family:sans-serif}
.seg:hover .tooltiptext{display:block}
.sg-struct{margin:.55rem 0 .15rem;border:1px solid #e7ecf3;border-radius:10px;padding:.5rem .6rem;background:#fbfdff}
.sg-head{font-size:.82rem;font-weight:700;color:#334155;margin:.1rem 0 .4rem;text-transform:uppercase;letter-spacing:.03em}
.sg-toolbar{display:flex;gap:.6rem;flex-wrap:wrap;align-items:center;margin:.2rem 0 .5rem}
.sg-dl{display:inline-block;border:1px solid #2563eb;background:#eff6ff;color:#1d4ed8;padding:.3rem .7rem;border-radius:6px;text-decoration:none;font-size:.8rem;font-weight:600}
.sg-dl:hover{background:#dbeafe}
.sg-map-svg{width:100%;overflow-x:auto;border:1px solid #e7ecf3;border-radius:8px;background:#fff;padding:4px 0}
.sg-map-svg svg{min-width:640px}
.feat-band{cursor:pointer}
.feat-band:hover{filter:brightness(1.06);stroke:#0f172a!important;stroke-width:1.6}
.sg-seq{font-family:ui-monospace,Consolas,monospace;font-size:12px;line-height:1.55;white-space:pre-wrap;word-break:break-all;max-height:230px;overflow:auto;border:1px solid #e7ecf3;border-radius:8px;padding:8px 10px;background:#fff;color:#334155;margin-top:.4rem}
.sg-seq b{font-weight:600}
.sg-table{width:100%;border-collapse:collapse;font-size:.78rem;margin-top:.5rem}
.sg-table th,.sg-table td{border-bottom:1px solid #eef2f7;padding:.28rem .45rem;text-align:left;vertical-align:top}
.sg-table th{color:#64748b;font-weight:600;background:#f8fafc}
.sg-role{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:middle}
.pill{display:inline-block;padding:1px 7px;border-radius:999px;font-size:.7rem;font-weight:600}
.pill-R1{background:#dcfce7;color:#166534}.pill-R2{background:#ede9fe;color:#5b21b6}.pill-ADP{background:#d1fae5;color:#065f46}.pill-CAP{background:#cffafe;color:#155e75}.pill-MASK{background:#e2e8f0;color:#334155}.pill-INL{background:#fef9c3;color:#713f12}.pill-POLY{background:#ffe4e6;color:#9f1239}.pill-BACK{background:#ede9fe;color:#5b21b6}.pill-INS{background:#ffedd5;color:#9a3412}
"""

_JS = """
(function(){
  var D = JSON.parse(document.getElementById('adapdata').textContent);
  var host = document.getElementById('adapp');
  var box = document.getElementById('adsearch');
  var countEl = document.getElementById('adcount');

  function esc(s){ return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function title(t){
    var label = t.nm || (D.roles[t.k] || t.k);
    if (t.k === 'insert'){ var m={'+':'sense','-':'antisense',':':'unstranded'}; return label + ' — ' + t.s + ' (' + (m[t.s]||t.s) + ')'; }
    var len = (t.s || '').length;
    var txt = label + ' (len ' + len + ')';
    if (t.n){ txt += ' -> rename tag {' + t.n + '}'; }
    return txt;
  }
  function seg(t, scheme){
    var cls = 'seg seg-' + t.k;
    var tip = title(t);
    // Copy the WHOLE scheme (not just this token) so a click anywhere on the
    // card yields the runnable scheme string, matching the "copy scheme" button.
    var copy = (scheme || t.s);
    return '<span class="'+cls+'" tabindex="0" data-copy="'+esc(copy)+'" data-tip="'+esc(tip)+'">'+esc(t.s)+'<span class="tooltiptext">'+esc(tip)+'</span></span>';
  }
  var ROLE_COLOR = {adp:'#A8E6CF', capture:'#B2EBF2', mask:'#e2e8f0', inline:'#FFD700', polytail:'#FFB3BA', back:'#C3B1E1', insert:'#FF8A65'};
  var ROLE_LBL = {adp:'adapter', capture:'UMI', mask:'mask', inline:'barcode', polytail:'poly-tail', back:'read-through', insert:'insert'};

  function mapSvg(c){
    var k = c.construct, size = k.size, feats = k.features;
    var w = 920, h = 66, x0 = 8, x1 = w - 8;
    function sx(p){ return x0 + (x1 - x0) * (p - 1) / Math.max(1, size - 1); }
    var s = '<div class="sg-map-svg"><svg viewBox="0 0 '+w+' '+h+'" width="100%">';
    s += '<line x1="'+x0+'" y1="34" x2="'+x1+'" y2="34" stroke="#cbd5e1" stroke-width="14" stroke-linecap="round"/>';
    s += '<text x="'+x0+'" y="14" font-size="10" fill="#94a3b8">5&prime;</text>';
    s += '<text x="'+x1+'" y="14" font-size="10" fill="#94a3b8" text-anchor="end">3&prime;</text>';
    for (var i=0;i<feats.length;i++){
      var f = feats[i];
      var a = sx(f.start), b = sx(f.end);
      if (b - a < 6){ var mid=(a+b)/2; a=mid-3; b=mid+3; }
      var col = ROLE_COLOR[f.role] || '#e0e0e0';
      s += '<g class="feat-band">';
      s += '<rect x="'+a+'" y="27" width="'+(b-a)+'" height="14" rx="3" fill="'+col+'" stroke="#94a3b8" stroke-width="1"/>';
      if (f.strand === '-'){ s += '<path d="M '+(a+7)+' 27 L '+(a-1)+' 34 L '+(a+7)+' 41 z" fill="'+col+'" stroke="#64748b"/>'; }
      else { s += '<path d="M '+(b-7)+' 27 L '+(b+1)+' 34 L '+(b-7)+' 41 z" fill="'+col+'" stroke="#64748b"/>'; }
      s += '</g>';
    }
    s += '<text x="'+(x0+x1)/2+'" y="60" font-size="10" fill="#64748b" text-anchor="middle">'+size+' bp &middot; linear</text>';
    s += '</svg></div>';
    s += '<div class="sg-leg" style="display:none"></div>';
    return s;
  }

  function seqView(c){
    var k = c.construct, size = k.size, feats = k.features;
    var colors = []; for (var i=0;i<size;i++) colors[i] = '#f8fafc';
    for (var i=0;i<feats.length;i++){ var f=feats[i]; var col=ROLE_COLOR[f.role]||'#e2e8f0'; for (var p=f.start-1;p<f.end;p++) colors[p]=col; }
    var s = '<div class="sg-seq">';
    for (var p=0;p<size;p++){ s += '<b style="background:'+colors[p]+'">'+k.seq[p]+'</b>'; if((p+1)%10===0) s+=' '; if((p+1)%50===0) s+='&#10;'; }
    s += '</div>';
    return s;
  }

  function pill(role, side){
    var p = '<span class="pill pill-'+role.toUpperCase()+'" style="margin-left:4px">'+ROLE_LBL[role]+'</span>';
    if (side) p += '<span class="pill pill-'+side+'" style="margin-left:4px">'+side+'</span>';
    return p;
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
    var s = '<div class="sg-struct">';
    s += '<div class="sg-toolbar"><span class="sg-head">Structure</span>';
    s += '<span style="margin-left:auto;display:flex;gap:.5rem">';
    s += '<a class="sg-dl" href="constructs/'+encodeURIComponent(c.name)+'.gb" download>Download .gb (GenBank)</a>';
    s += '</span></div>';
    s += mapSvg(c);
    s += seqView(c);
    s += featTable(c);
    s += '</div>';
    return s;
  }

  function card(c){
    var alias = c.alias ? ' data-alias="1"' : '';
    var cmd = 'cutseq -A ' + c.name;
    var points = '';
    if (c.points && c.points.length){ points = '<ul class="card-points">' + c.points.map(function(p){return '<li>'+esc(p)+'</li>';}).join('') + '</ul>'; }
    var s = '';
    s += '<details class="adapter-card"'+alias+'>';
    s += '<summary data-copy="'+esc(c.scheme)+'"><span class="card-name">'+esc(c.name)+'</span><span class="card-desc">'+esc(c.desc)+'</span>'+(c.alias?'<span class="card-desc">[alias of '+esc(c.alias)+']</span>':'')+'</summary>';
    s += '<div class="card-body"><div class="seg-row">' + c.tokens.map(function(t){ return seg(t, c.scheme); }).join('') + '</div>';
    s += structBlock(c);
    s += '<div class="card-actions">';
    s += '<button class="copy-btn" data-copy="'+esc(c.scheme)+'">copy scheme</button>';
    s += '<button class="copy-btn" data-copy="'+esc(cmd)+'">copy: '+esc(cmd)+'</button>';
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

  document.addEventListener('click', function(e){
    var el = e.target.closest ? e.target.closest('[data-copy]') : null;
    if (!el){ return; }
    var val = el.getAttribute('data-copy');
    function done(){
      var old = el.getAttribute('data-label');
      if (!old){ el.setAttribute('data-label', el.textContent); el.textContent = 'copied!'; setTimeout(function(){ el.textContent = el.getAttribute('data-label'); el.removeAttribute('data-label'); }, 900); }
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
  .ad-toolbar{background:#1e293b;border-color:#334155}
  .ad-search{background:#0f172a;color:#e2e8f0;border-color:#334155}
  .ad-count{color:#94a3b8}
  .ad-legend{background:#0f172a;border-color:#334155}
  .adleg{color:#cbd5e1}
  .adapter-card{background:#0f172a;border-color:#334155}
  .card-name{color:#e2e8f0}
  .card-desc{color:#94a3b8}
  .card-points{color:#cbd5e1}
  .sg-struct{background:#0f172a;border-color:#334155}
  .sg-head{color:#cbd5e1}
  .sg-map-svg{background:#0f172a;border-color:#334155}
  .sg-seq{background:#0f172a;color:#e2e8f0;border-color:#334155}
  .sg-table th{background:#1e293b;color:#cbd5e1}
  .sg-table th,.sg-table td{border-color:#334155}
  .copy-btn{background:#1e293b;color:#e2e8f0;border-color:#334155}
  .copy-btn:hover{background:#1e293b;border-color:#60a5fa;color:#e2e8f0}
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
    "    <div class=\"search-box\">\n"
    "      <input id=\"site-search\" class=\"site-search-input\" type=\"search\" placeholder=\"Search docs&hellip;\" autocomplete=\"off\">\n"
    "      <div id=\"site-search-results\" class=\"site-search-results\" role=\"listbox\" hidden></div>\n"
    "    </div>\n"
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