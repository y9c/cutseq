#!/usr/bin/env python3
"""Post-migration gate for the cutseq docs static-site.

Independent, read-only verification of the generated static site under ``docs/``.
The site is now fully static HTML (the ``docs/*.md`` sources have been removed),
so content fidelity is checked **against the generated HTML itself**: distinctive
mid-body technical phrases that must survive regeneration are asserted directly
in each page's visible text.  No source/page files are touched; this module never
edits anything and never needs a browser.

Every assertion is derived from the spec guarantees:
  * canonical prose phrases appear verbatim in each generated page
  * markdown tables / callouts (``{: .note}`` / ``{: .warning}``) / images /
    fenced code blocks map to ``<table>`` / ``.callout`` / ``<img>`` / ``<pre><code>``
  * cross-page relative links resolve to real files (and never ``.md``)
  * a shared shell (``<nav>`` with the 4 page links, GA id G-EFHY8WFMHK,
    footer "Copyright © 2021-2026, Chang Ye", a "View on GitHub" aux link)
  * a client search index (``assets/js/search.json``) parses with >=4 entries
  * ``adapters.html`` keeps a working viewer (its ``adapdata`` JSON parses and has
    ``schemes``; its viewer ``<script>`` is at least statically sane)
  * ``index.html`` exists at site root and is a complete HTML document
  * no leftover ``.md`` hrefs / remnants anywhere in the generated HTML

Nothing here depends on implementation-specific class names or element ids beyond
``.callout``, ``nav`` and ``search.json``.
"""

from __future__ import annotations

import json
import re
from html import unescape as _unescape
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"

# The four prose pages carry the content we verify.  Ground truth is the emitted
# HTML itself (the markdown sources no longer exist).  The former "spatial" page
# was folded into adapters.html (see the DBiT-seq section), so it is no longer a
# standalone page.
PROSE_PAGES = ["index", "quickstart", "rename"]
# Every page that must share the shell (including the standalone adapters viewer).
ALL_PAGES = PROSE_PAGES + ["adapters"]
# Cross-page targets the shared shell / content must be able to reach.
CROSS_PAGE = ["index.html", "quickstart.html", "rename.html", "adapters.html"]

# Canonical substrings picked from the generated HTML (mid-body, technical,
# stable, pure plain text -- no markdown emphasis/backtick markers, no HTML
# entities).  Each must survive inside the generated page's visible text.
CANON = {
    "index": [
        "Automatically cut adapter / barcode / UMI from NGS data, in the right order.",
        "NGS library prep can need many sequential trimming and extraction steps.",
        "so the whole pipeline is one parameter and one command.",
        "capture and rename inline barcodes / UMIs into read names.",
    ],
    "quickstart": [
        "The scheme string describes the structure of your library.",
        ( "A scheme is compiled once into a lazy graph of native cutadapt "
          "modifiers and executed in a single pass." ),
        "Defaults reproduce legacy output byte-for-byte.",
        "fixed uppercase run directly adjacent to a lowercase edge",
    ],
    "rename": [
        "By default, CutSeq appends every captured UMI / inline barcode to the read name",
        "The template is a mix of literal text",
        "Any capture may be wrapped in functions.",
        "all produce a clear error message at build time",
    ],
}

# Minimum structural marks that must be present in each page's main region,
# derived from the current generated HTML (so a future regeneration that DROPS a
# table / callout / image / code block is flagged).
EXPECTED_STATS = {
    "index": {"tables": 1, "callouts": 1, "images": 1, "pre": 0},
    "quickstart": {"tables": 0, "callouts": 0, "images": 1, "pre": 6},
    "rename": {"tables": 2, "callouts": 2, "images": 0, "pre": 3},
}

GA_ID = "G-EFHY8WFMHK"
FOOTER_YEAR = "2021-2026"
FOOTER_OWNER = "Chang Ye"
GITHUB_TEXT = "View on GitHub"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _norm(text: str) -> str:
    """Collapse all runs of whitespace to single spaces (handles md soft-wraps)."""
    return re.sub(r"\s+", " ", text).strip()


def _read(rel: Path) -> str:
    return (DOCS / rel).read_text(encoding="utf-8")


def _extract_main(h: str) -> str:
    m = re.search(r"<main\b[^>]*>(.*?)</main>", h, re.S)
    if m:
        return m.group(1)
    m = re.search(r"<body[^>]*>(.*)</body>", h, re.S)
    return m.group(1) if m else h


def _page_text(h: str) -> str:
    """Visible text of the main content region (tags stripped, entities unescaped)."""
    frag = _extract_main(h)
    frag = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", " ", frag, flags=re.S | re.I)
    frag = re.sub(r"<[^>]+>", " ", frag)
    return _norm(_unescape(frag))


def _callout_count(frag: str) -> int:
    """Count elements with an exact ``callout`` class token (excludes .callout-title)."""
    n = 0
    for cls in re.findall(r'class="([^"]*)"', frag):
        if re.search(r"(?<![A-Za-z0-9-])callout(?![A-Za-z0-9-])", cls):
            n += 1
    return n


def _main_stats(h: str) -> dict:
    frag = _extract_main(h)
    return {
        "tables": len(re.findall(r"<table\b", frag)),
        "callouts": _callout_count(frag),
        "images": len(re.findall(r"<img\b", frag)),
        "pre": len(re.findall(r"<pre\b", frag)),
    }


def _relative_refs(h: str) -> list:
    """Static relative ``href``/``src`` targets (skips JS-built and absolute)."""
    refs = []
    for r in re.findall(r'(?:href|src)="([^"]+)"', h):
        if r.startswith(("http://", "https://", "//", "#", "mailto:", "tel:", "/", "data:")):
            continue
        if "'" in r or "(" in r or "+" in r:  # JS expression, not a static link
            continue
        target = r.split("#")[0].split("?")[0]
        if target:
            refs.append(target)
    return refs


def _nav_block(h: str) -> str:
    m = re.search(r"<nav\b[^>]*>(.*?)</nav>", h, re.S)
    return m.group(1) if m else ""


def _viewer_js(h: str) -> str:
    """The adapters viewer JS (the bare <script> that reads the adapdata JSON)."""
    for blk in re.findall(r"<script>(.*?)</script>", h, re.S):
        if "JSON.parse" in blk and "adapdata" in blk:
            return blk
    return ""


def _strip_js(js: str) -> str:
    s = re.sub(r"/\*.*?\*/", "", js, flags=re.S)
    s = re.sub(r"//[^\n]*", "", s)
    s = re.sub(r"'(?:\\.|[^'\\])*'", "''", s)
    s = re.sub(r'"(?:\\.|[^"\\])*"', '""', s)
    s = re.sub(r"`(?:\\.|[^`\\])*`", "``", s)
    return s


# --------------------------------------------------------------------------- #
# 1. content fidelity (prose pages -- ground truth is the generated HTML)
# --------------------------------------------------------------------------- #
def test_content_fidelity():
    """Every canonical phrase survives into the generated page's visible text."""
    for page in PROSE_PAGES:
        text = _page_text(_read(Path(page + ".html")))
        for sub in CANON[page]:
            assert sub in text, f"[{page}] content lost in {page}.html: {sub!r}"


# --------------------------------------------------------------------------- #
# 2. structural mark conversion (tables / callouts / images / code)
# --------------------------------------------------------------------------- #
def test_tables_callouts_images_code_converted():
    """Markdown tables -> <table>, callouts -> .callout, images -> <img>,
    fenced code -> <pre><code>.  Ground truth is the expected structural counts
    derived from the current generated HTML."""
    for page in PROSE_PAGES:
        stats = _main_stats(_read(Path(page + ".html")))
        exp = EXPECTED_STATS[page]
        name = page + ".html"
        assert stats["tables"] >= exp["tables"], f"{name}: dropped tables"
        if exp["tables"]:
            assert stats["tables"] > 0, f"{name}: expected tables but none present"
        assert stats["callouts"] >= exp["callouts"], f"{name}: dropped callouts"
        if exp["callouts"]:
            assert stats["callouts"] > 0, f"{name}: expected callouts but none present"
        assert stats["images"] >= exp["images"], f"{name}: dropped images"
        if exp["images"]:
            assert stats["images"] > 0, f"{name}: expected images but no <img>"
        assert stats["pre"] >= exp["pre"], f"{name}: dropped code blocks"
        if exp["pre"]:
            assert "<pre><code" in _extract_main(_read(Path(page + ".html"))), \
                f"{name}: code not wrapped in <pre><code>"


# --------------------------------------------------------------------------- #
# 3. cross-page relative links resolve
# --------------------------------------------------------------------------- #
def test_relative_links_resolve():
    """All static relative links on every page resolve to an existing file."""
    for page in ALL_PAGES:
        h = _read(Path(page + ".html"))
        for target in _relative_refs(h):
            assert not target.endswith(".md"), f"{page}.html links to .md: {target!r}"
            assert (DOCS / target).is_file(), f"{page}.html broken link: {target!r}"
    # the specific cross-page targets must all exist at site root
    for target in CROSS_PAGE:
        assert (DOCS / target).is_file(), f"missing cross-page target: {target}"


# --------------------------------------------------------------------------- #
# 4. shared shell on every page
# --------------------------------------------------------------------------- #
def test_shared_shell():
    """Every page has the shared nav (4 page links, no spatial), GA, footer, GitHub aux link."""
    for page in ALL_PAGES:
        h = _read(Path(page + ".html"))
        nav = _nav_block(h)
        assert nav, f"{page}.html missing <nav>"
        # exactly the 4 current nav items; the old "Spatial (DBiT-seq)" item is gone
        items = re.findall(r'<li class="nav-item[^"]*">', nav)
        assert len(items) == 4, f"{page}.html nav has {len(items)} items, expected 4"
        for target in CROSS_PAGE:
            assert f'href="{target}"' in nav, f"{page}.html nav missing link to {target}"
        assert "spatial.html" not in nav, f"{page}.html nav still links to removed spatial.html"
        assert GA_ID in h, f"{page}.html missing GA id {GA_ID}"
        assert FOOTER_YEAR in h and FOOTER_OWNER in h, f"{page}.html missing footer"
        assert GITHUB_TEXT in h, f"{page}.html missing 'View on GitHub' aux link"


# --------------------------------------------------------------------------- #
# 5. search index
# --------------------------------------------------------------------------- #
def test_search_index():
    """assets/js/search.json parses to >=4 entries, each with title+url+text.

    The former "Spatial (DBiT-seq) libraries" entry is gone; its DBiT keywords are
    folded into the (single) Adapter Schemes entry, which points at the DBiT section
    anchor on adapters.html."""
    p = DOCS / "assets" / "js" / "search.json"
    assert p.is_file(), f"missing {p}"
    data = json.loads(p.read_text(encoding="utf-8"))
    assert isinstance(data, list), "search.json must be a JSON array"
    assert len(data) >= 4, f"search.json has {len(data)} entries, need >=4"
    urls = []
    for i, entry in enumerate(data):
        assert isinstance(entry, dict), f"entry {i} not an object"
        for field in ("title", "url", "text"):
            v = entry.get(field)
            assert isinstance(v, str) and v.strip(), f"entry {i} missing non-empty {field}"
        assert not entry["url"].endswith(".md"), f"entry {i} url still points to .md"
        urls.append(entry["url"])
    # the removed spatial page must not be a search target, but its DBiT content must
    # still be findable via the Adapter Schemes entry (now inside the DBiT section).
    assert "spatial.html" not in urls, "search index still points at removed spatial.html"
    adapters = next(e for e in data if e["title"] == "Adapter Schemes")
    assert adapters["url"].startswith("adapters.html"), \
        "Adapter Schemes entry should target adapters.html"
    folded = adapters["text"]
    for kw in ("DBiT-seq", "BarcodeB", "41-nt peak", "barcode arm"):
        assert kw.lower() in folded.lower(), f"Adapter Schemes search text lost DBiT keyword {kw!r}"


# --------------------------------------------------------------------------- #
# 6. adapters viewer
# --------------------------------------------------------------------------- #
def test_adapters_viewer_json():
    """adapters.html still embeds a parseable adapdata JSON with 'schemes'."""
    h = _read(Path("adapters.html"))
    m = re.search(r'<script type="application/json"[^>]*id="adapdata"[^>]*>(.*?)</script>',
                  h, re.S)
    assert m, "adapters.html missing adapdata JSON block"
    data = json.loads(m.group(1))
    assert "schemes" in data, "adapdata JSON missing 'schemes'"
    assert len(data["schemes"]) >= 1, "adapdata has no schemes"


def test_adapters_viewer_js_sanity():
    """The adapters viewer <script> is non-trivial and statically balanced."""
    js = _viewer_js(_read(Path("adapters.html")))
    assert len(js) > 200, "viewer <script> missing/too short"
    assert "adapdata" in js, "viewer must read the adapdata payload"
    assert "JSON.parse" in js, "viewer must parse adapdata"
    stripped = _strip_js(js)
    assert stripped.count("{") == stripped.count("}"), "viewer JS: unbalanced braces"
    assert stripped.count("(") == stripped.count(")"), "viewer JS: unbalanced parens"
    assert stripped.count("[") == stripped.count("]"), "viewer JS: unbalanced brackets"


# --------------------------------------------------------------------------- #
# 6b. folded spatial (DBiT-seq) section inside adapters.html
# --------------------------------------------------------------------------- #
def test_adapters_dbit_folded():
    """The former spatial.html content is folded into adapters.html as a detailed
    section (after the discard section), and spatial.html no longer exists."""
    h = _read(Path("adapters.html"))
    main = _extract_main(h)
    text = _page_text(h)
    # distinctive DBiT prose that must survive the fold
    for sub in [
        "DBiT-seq puts a row of fixed-length oligos",
        "the fixed barcode arm occupies 110 nt",
        "R1 lengths reflect the real insert-length spread",
        "sequences start downstream of the priming site.",
    ]:
        assert sub in text, f"adapters.html lost folded DBiT phrase: {sub!r}"
    # the section anchor + the deep per-paragraph anchors survive
    assert '<section id="spatial-dbit-seq">' in main, "adapters.html missing DBiT section wrapper"
    assert 'id="5-interpreting-the-read-length-distribution"' in main, \
        "adapters.html lost the R2 41-nt peak anchor"
    # the DBITSEQ-scheme/token table (Read/What-it-sequences + Token/Side/Meaning) still present
    assert main.count("<table>") >= 2, "adapters.html lost the DBiT tables"
    assert "DBITSEQ" in main, "adapters.html lost the DBITSEQ scheme reference"
    # the old standalone page is gone and nothing links to it
    assert not (DOCS / "spatial.html").exists(), "spatial.html still exists"
    for page in ALL_PAGES:
        assert "spatial.html" not in _read(Path(page + ".html")), \
            f"{page}.html still references removed spatial.html"


# --------------------------------------------------------------------------- #
# 7. site root / index
# --------------------------------------------------------------------------- #
def test_index_at_site_root():
    """index.html is the site root document and is a complete HTML page."""
    p = DOCS / "index.html"
    assert p.is_file(), "docs/index.html does not exist"
    h = p.read_text(encoding="utf-8")
    assert re.search(r"^\s*<!doctype html>\s*$", h, re.I | re.M), "index.html missing <!doctype html>"
    assert re.search(r"</html>\s*$", h), "index.html missing closing </html>"
    assert "<title>" in h and "<body" in h, "index.html not a complete page"


# --------------------------------------------------------------------------- #
# 8. no .md link remnants
# --------------------------------------------------------------------------- #
def test_no_md_link_remnants():
    """No generated page may reference markdown files from links."""
    leaked = []
    for page in ALL_PAGES:
        h = _read(Path(page + ".html"))
        if re.search(r'(?:href|src)="[^"]*\.md"', h):
            leaked.append(page + ".html")
        for t in ("adapters.md", "quickstart.md", "rename.md", "spatial.md", "index.md"):
            if re.search(rf'href="[^"]*{re.escape(t)}"', h):
                leaked.append(f"{page}.html ({t})")
    assert not leaked, f"markdown link remnants: {leaked}"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
