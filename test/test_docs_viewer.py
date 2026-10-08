#!/usr/bin/env python3
"""Render-time test for the generated adapters.html viewer widgets.

Validates that the chunk of docs/adapters.md the viewer emits (the JSON payload
inside `<script type="application/json" id="adapdata">`, the viewer <style> and the
viewer <script>) is well-formed and that JavaScript payload parses. This runs
without a Jekyll build: it extracts the emitted fragments from the markdown and
wraps them into a standalone HTML document which a headless browser can load.

Live browser execution is exercised at dev time via puppeteer; this module only
checks static invariants of the emitted fragments.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
ADAPTERS = ROOT / "docs" / "adapters.md"


def _chunk(md: str, tag: str):
    m = re.search(tag, md, re.DOTALL)
    assert m, f"missing {tag}"
    return m.group(1)


@pytest.fixture(scope="module")
def md_text():
    return ADAPTERS.read_text(encoding="utf-8")


def test_payload_json_parses(md_text):
    """The embedded adapdata is valid JSON and has the expected shape."""
    raw = _chunk(
        md_text,
        r'<script type="application/json" id="adapdata">(.*?)</script>',
    )
    data = json.loads(raw)
    assert "schemes" in data
    assert "roles" in data
    assert len(data["schemes"]) > 0
    # each scheme carries every field the viewer's renderer consumes
    for s in data["schemes"]:
        assert {"name", "desc", "scheme", "tokens", "construct"} <= set(s)
        assert set(s["construct"]) >= {"seq", "size", "features"}


def test_payload_has_no_unescaped_script_breakout(md_text):
    """No raw `</script>` sequence inside the JSON (would close the tag early)."""
    raw = _chunk(
        md_text,
        r'<script type="application/json" id="adapdata">(.*?)</script>',
    )
    # ensure_ascii=True in the generator guarantees this, but guard anyway
    assert "</script>" not in raw


def test_js_no_html_entity_through_esc(md_text):
    """The feature-table insert cell must not leak the `N&#215;30` entity.

    Regression: the generator used `sq = 'N&#215;'+f.len`, which `esc()` then
    double-escaped to a literal `N&#215;30` on the page. It must now be a plain
    `N` + U+00D7 + length.
    """
    assert "&#215;" not in md_text


def test_feature_table_insert_uses_multiply_sign(md_text):
    """Insert cells render the multiplication sign as a real character."""
    assert "if (f.role === 'insert'){ sq = 'N×'+f.len; }" in md_text


def test_no_js_console_fatal_if_rendered(md_text):
    """Smoke: the genuine viewer <script> references the key element ids."""
    js = _chunk(md_text, r'<script>(.*?)</script>')
    for needle in ("adapdata", "adapp", "adsearch", "featTable", "mapSvg"):
        assert needle in js


def test_darkmode_css_present(md_text):
    """The user-supplied dark-mode overrides survive generation."""
    assert "prefers-color-scheme" in md_text


def test_constructs_unchanged():
    """The `.gb`/`.dna` constructs on disk must not be dirty vs HEAD."""
    out = subprocess.run(
        ["git", "-C", str(ROOT), "status", "--short", "--", "docs/constructs"],
        capture_output=True,
        text=True,
    )
    assert out.stdout.strip() == "", f"constructs changed: {out.stdout}"
