#!/usr/bin/env python3
"""Generate docs/adapters.md as an interactive adapter viewer.

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


def _token_json(tok, ordinal):
    """A token -> {k, s, n} for the viewer (k=kind, s=display text, n=capture ord)."""
    kind = tok.kind
    if kind == "capture":
        return {"k": kind, "s": _display_len(kind, tok.value), "n": ordinal}
    if kind == "mask":
        return {"k": kind, "s": _display_len(kind, tok.value), "n": None}
    if kind == "insert":
        return {"k": kind, "s": tok.value, "n": None}
    if kind == "polytail":
        return {"k": kind, "s": f"{tok.value}...{tok.value}", "n": None}
    return {"k": kind, "s": tok.value, "n": None}


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


def _build_construct(tokens):
    """Assemble the top-strand molecule 5'->3' = R1-side + insert + R2-side.

    The R2-side of a scheme is already written top-strand 5'->3', so appending
    tokens in order (with an insert placeholder at the marker) yields the full
    construct. Each token becomes an annotated feature with 1-based coords.
    """
    side = "R1"
    parts = []
    feats = []
    cursor = 0
    for tok in tokens:
        kind = tok.kind
        if kind == "insert":
            length = INSERT_PH
            text = "N" * length
            name = "library insert (placeholder)"
            ftype = "misc_feature"
            strand = "+"
        elif kind == "capture":
            length = int(tok.value)
            text = "N" * length
            name = f"UMI/capture N{length} (renamed to {tok.label or ('{' + str(tok.index) + '}')})" if tok.index else f"UMI/capture N{length}"
            ftype = "misc_feature"
            strand = "+"
        elif kind == "mask":
            length = int(tok.value)
            text = "N" * length
            name = f"mask X{length} (trimmed, not kept)"
            ftype = "misc_feature"
            strand = "+"
        elif kind == "adp":
            text = tok.value
            length = len(text)
            name = "adapter / primer"
            ftype = "misc_feature"
            strand = "+"
        elif kind == "inline":
            text = tok.value
            length = len(text)
            name = "inline barcode"
            ftype = "misc_feature"
            strand = "+"
        elif kind == "polytail":
            length = POLY_RUN
            text = tok.value * length
            name = f"poly-{tok.value} tail"
            ftype = "misc_feature"
            strand = "+"
        elif kind == "back":
            text = tok.value
            length = len(text)
            name = "3' read-through adapter"
            ftype = "misc_feature"
            strand = "-"
        else:
            continue
        start = cursor + 1
        end = cursor + length
        feats.append({"name": name, "type": ftype, "role": kind, "strand": strand,
                      "side": side, "start": start, "end": end, "len": length,
                      "seq": text})
        parts.append(text)
        cursor += length
        # switch the read side snapshot once we pass the insert marker
        if kind == "insert":
            side = "R2"
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
        for tok in raw_tokens:
            tok_ord = None
            if tok.kind in ("capture", "inline"):
                ordinal += 1
                tok_ord = ordinal
            tokens.append(_token_json(tok, tok_ord))
        # annotate capture ordinals back onto the raw tokens for the construct
        idx = 0
        for tok in raw_tokens:
            if tok.kind in ("capture", "inline"):
                idx += 1
                tok.index = idx
        construct = _build_construct(raw_tokens)
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
    var label = D.roles[t.k] || t.k;
    if (t.k === 'insert'){ var m={'+':'sense','-':'antisense',':':'unstranded'}; return 'insert marker ' + t.s + ' (' + (m[t.s]||t.s) + ')'; }
    var len = (t.s || '').length;
    var txt = label + ' (len ' + len + ')';
    if (t.n){ txt += ' -> rename tag {' + t.n + '}'; }
    return txt;
  }
  function seg(t){
    var cls = 'seg seg-' + t.k;
    var tip = title(t);
    return '<span class="'+cls+'" tabindex="0" data-copy="'+esc(t.s)+'" data-tip="'+esc(tip)+'">'+esc(t.s)+'<span class="tooltiptext">'+esc(tip)+'</span></span>';
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
      if (f.role === 'insert'){ sq = 'N&#215;'+f.len; }
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
    s += '<a class="sg-dl" href="constructs/'+encodeURIComponent(c.name)+'.dna" download>Download .dna (SnapGene)</a>';
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
    s += '<summary><span class="card-name">'+esc(c.name)+'</span><span class="card-desc">'+esc(c.desc)+'</span>'+(c.alias?'<span class="card-desc">[alias of '+esc(c.alias)+']</span>':'')+'</summary>';
    s += '<div class="card-body"><div class="seg-row">' + c.tokens.map(seg).join('') + '</div>';
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
    if (navigator.clipboard && navigator.clipboard.writeText){ navigator.clipboard.writeText(val).then(done); }
    else { var t = document.createElement('textarea'); t.value = val; document.body.appendChild(t); t.select(); try{ document.execCommand('copy'); }catch(err){}; document.body.removeChild(t); done(); }
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
    from datetime import date
    return date.today().strftime("%d-%b-%Y").upper()


def _write_exports(name, construct, constructs_dir):
    """Write a GenBank (.gb) and a SnapGene (.dna) file for one scheme."""
    constructs_dir.mkdir(parents=True, exist_ok=True)
    gb = _genbank(name, "", "", construct["seq"], construct["features"])
    (constructs_dir / f"{name}.gb").write_text(gb, encoding="utf-8")
    try:
        from sgffp import SgffObject, SgffWriter
    except Exception as e:  # pragma: no cover
        print(f"  [warn] sgffp unavailable, skipping .dna for {name}: {e}")
        return
    obj = SgffObject.new(construct["seq"], topology="linear")
    for f in construct["features"]:
        try:
            obj.add_feature(f["name"], f["type"], f["start"] - 1, f["end"])
        except Exception as e:
            print(f"  [warn] feature skipped for {name}: {e}")
    SgffWriter.to_file(obj, str(constructs_dir / f"{name}.dna"))


def main():
    adapters_toml = project_root / "cutseq" / "adapters.toml"
    out_md = project_root / "docs" / "adapters.md"
    schemes = _load_schemes()
    if not schemes:
        sys.exit("Error: no schemes loaded from adapters.toml")

    payload = {"schemes": schemes, "roles": {k: v[0] for k, v in ROLES.items()},
               "legend": _legend_html()}
    data_json = json.dumps(payload, ensure_ascii=True)

    # write downloadable construct files (.gb + .dna)
    constructs_dir = project_root / "docs" / "constructs"
    for s in schemes:
        _write_exports(s["name"], s["construct"], constructs_dir)
    print(f"wrote constructs for {len(schemes)} schemes to {constructs_dir}")

    payload = {"schemes": schemes, "roles": {k: v[0] for k, v in ROLES.items()},
               "legend": _legend_html()}
    data_json = json.dumps(payload, ensure_ascii=True)

    md = []
    md.append("---\ntitle: Adapter Schemes\nnav_order: 3\n---\n\n")
    md.append("# Adapter Schemes\n\n")
    md.append("CutSeq bundles built-in library schemes — list them in a terminal "
              "with `cutseq --list-adapters`, and use a name (or a custom grammar "
              "string) via `-A/--adapter-scheme` (`-a` is an alias).\n\n")
    md.append("This page is a live viewer: **search** to filter, **hover** a colored "
              "segment to see its role, **click** any segment or button to copy it.\n\n")
    md.append(_legend_html())
    md.append('\n<div class="ad-toolbar">')
    md.append('<input id="adsearch" class="ad-search" type="search" '
              'placeholder="Search schemes…" />')
    md.append(f'<span id="adcount" class="ad-count">{len(schemes)} schemes</span>')
    md.append('</div>\n')
    md.append('<div id="adlegend-host" style="display:none"></div>')
    md.append('<div id="adapp"></div>\n')
    md.append(f'<script type="application/json" id="adapdata">{data_json}</script>\n')
    md.append('<style>' + _CSS + '</style>\n')
    md.append('<script>' + _JS + '</script>\n')

    md.append("\n## Inline-barcode auto-detection\n\n")
    md.append("Inline barcodes are written in lowercase (`acgt...`). If you write "
              "one in uppercase, it merges with the adjacent sequencing primer into "
              "one adapter run. CutSeq checks the two outermost adapters of a scheme "
              "against a curated database of Illumina / BGI (MGI) sequencing primers "
              "(`cutseq --list-primers`) and reclassifies any fixed uppercase "
              "sequence adjacent to (or between) recognized primers as an inline "
              "barcode, with a warning. Custom schemes without known primers are "
              "never altered. Disable with `--no-auto-inline`.\n\n")
    md.append("## Discarding reads with a reason tag\n\n")
    md.append("A single `-d` output collects every rejected read with the reason in "
              "the name (`reason=too_short`, `reason=too_many_n`, `reason=low_quality`, "
              "`reason=no_barcode`), gated by `-m/--min-length`, `--max-n`, "
              "`--min-avg-quality` and `--ensure-inline-barcode`. When no `-d`/`-O` is "
              "given, discard files are auto-named from the inputs.\n")

    out_md.write_text("".join(md), encoding="utf-8")
    print(f"wrote {out_md}: {len(schemes)} schemes")


if __name__ == "__main__":
    main()