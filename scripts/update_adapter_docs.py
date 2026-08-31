#!/usr/bin/env python3
import sys
from pathlib import Path

# Add project root to sys.path to allow importing cutseq.grammar
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

try:
    from cutseq.grammar import tokenize
except ImportError as e:
    print(f"Error importing tokenize from cutseq.grammar: {e}")
    print(f"Current sys.path: {sys.path}")
    print(f"Project root: {project_root}")
    sys.exit(1)

try:
    import tomllib  # Python 3.11+
except ImportError:
    try:
        import tomli as tomllib  # Python < 3.11
    except ImportError:
        print(
            "Error: 'tomli' is not installed. Please install it with 'pip install tomli' for Python < 3.11."
        )
        sys.exit(1)

COLOR_PALETTE = {
    "adp": "#A8E6CF",
    "inline": "#FFD700",
    "capture": "#B2EBF2",
    "mask": "#DCDCDC",
    "polytail": "#FFB3BA",
    "insert_bg": "#FF6F61",  # A distinct background for the insert indicator
    "default_seq": "#E0E0E0",
}

# Token kinds that are rendered as the highlighted strand/insert hexagon.
_INSERT_KINDS = {"insert"}


def _token_color(kind):
    return COLOR_PALETTE.get(kind, COLOR_PALETTE["default_seq"])


def render_scheme_html(scheme_raw):
    """Render a grammar scheme as a color-coded, copyable HTML block."""
    tokens = tokenize(scheme_raw)
    raw_value = scheme_raw

    parts = [
        '<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;">',
        f'<div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: {raw_value}" data-scheme="{raw_value}">',
    ]

    for tok in tokens:
        if tok.kind in _INSERT_KINDS:
            # Render the R1|R2 split marker as a hexagon like the strand glyph.
            ch = tok.value
            parts.append(
                f'<div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: {COLOR_PALETTE["insert_bg"]}; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">{ch}</span></div>'
            )
        else:
            seq_val = _token_value(tok)
            color = _token_color(tok.kind)
            parts.append(
                f'<span style="background-color: {color}; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">{seq_val}</span>'
            )

    parts.append("</div>")  # Close flex container
    parts.append(
        '<div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div>'
    )
    parts.append("</div>\n")  # Close adapter-scheme div
    return "".join(parts)


def _token_value(tok):
    """Render a token's value for display (expanding N/X lengths to runs)."""
    kind = tok.kind
    v = tok.value
    if kind in ("capture", "mask"):
        return {"capture": "N", "mask": "X"}[kind] * v
    return str(v)


def main():
    adapters_toml_path = project_root / "cutseq" / "adapters.toml"
    adapters_md_path = project_root / "docs" / "adapters.md"

    print(f"Adapters TOML path: {adapters_toml_path}")
    print(f"Adapters MD path: {adapters_md_path}")

    if not adapters_toml_path.exists():
        print(
            f"Error: {adapters_toml_path} not found. Make sure the script is run from the project root or the path is correct."
        )
        sys.exit(1)

    try:
        with open(adapters_toml_path, "rb") as f:
            adapters_data = tomllib.load(f)
    except Exception as e:
        print(f"Error parsing {adapters_toml_path}: {e}")
        sys.exit(1)

    # --- Generate Markdown Content ---
    all_markdown_parts = []
    all_markdown_parts.append("---\ntitle: Adapter Schemes\nnav_order: 3\n---\n\n")
    all_markdown_parts.append("# Adapter Schemes\n\n")
    all_markdown_parts.append(
        "CutSeq supports a variety of built-in adapter schemes for common NGS library types. You can list all available schemes in your terminal with:\n\n"
    )
    all_markdown_parts.append("```bash\ncutseq --list-adapters\n```\n\n")
    all_markdown_parts.append(
        "Use the adapter name (or a custom grammar scheme string) with `-A/--adapter-scheme`. Note `-a` is accepted as an alias of `-A`.\n\n"
    )
    all_markdown_parts.append("## Example: Built-in Schemes\n\n")
    all_markdown_parts.append(
        "- **SMALLRNA**: Small RNA libraries, double ligation, forward orientation\n"
        "- **INLINE**: Custom barcoded libraries, dual UMI, inline barcode\n"
        "- **TAKARAV3**: SMARTer Stranded Total RNA-Seq Kit v3\n"
        "- **STRANDED**: Stranded RNA libraries\n\n"
    )
    all_markdown_parts.append(
        "See below for a comprehensive guide to each supported adapter pattern, including copyable scheme blocks and usage notes.\n\n---\n\n"
    )
    all_markdown_parts.append(
        "## Inline barcode auto-detection\n\n"
        "Inline barcodes are written in lowercase (`acgt...`). If you write one "
        "in uppercase, it merges with the adjacent sequencing primer into one "
        "adapter run. CutSeq checks the two outermost adapters of a scheme "
        "against a curated database of Illumina / BGI (MGI) sequencing primers "
        "(`cutseq --list-primers`) and reclassifies any fixed uppercase sequence "
        "adjacent to (or between) recognized primers as an inline barcode, with "
        "a warning. Custom schemes without known primers are never altered. "
        "Disable with `--no-auto-inline`.\n\n"
    )
    all_markdown_parts.append(
        "## Discarding reads with a reason tag\n\n"
        "A single discard output captures all reads that fail QC, with the reason "
        "stored in the read name (`reason=too_short`, `reason=too_many_n`, "
        "`reason=low_quality` or `reason=no_barcode`):\n\n"
        "```bash\n"
        "cutseq -A SMALLRNA -d discard_R1.fq.gz discard_R2.fq.gz \\\n"
        "    --min-length 20 --max-n 5 --min-avg-quality 20 \\\n"
        "    test_R1.fq.gz test_R2.fq.gz\n"
        "```\n\n"
        "- `reason=too_short` — shorter than `--min-length` after trimming\n"
        "- `reason=too_many_n` — exceeds `--max-n`\n"
        "- `reason=low_quality` — mean Phred quality below `--min-avg-quality`\n"
        "- `reason=no_barcode` — missing an expected inline barcode (only with "
        "`--ensure-inline-barcode`, which requires inline barcodes in the scheme)\n\n"
        "When no `-d`/`-O` is given, discard files are auto-named from the input "
        "files.\n\n---\n\n"
    )

    for adapter_key, adapter_info in adapters_data.items():
        if not isinstance(adapter_info, dict) or not all(
            k in adapter_info for k in ["scheme", "description_name", "points"]
        ):
            print(
                f"Warning: Skipping adapter '{adapter_key}' due to missing 'scheme', 'description_name', or 'points' fields."
            )
            continue

        all_markdown_parts.append(
            f"### {adapter_key} ({adapter_info['description_name']})\n\n"
        )

        scheme_raw = adapter_info["scheme"]
        all_markdown_parts.append(render_scheme_html(scheme_raw))

        # Bullet Points
        if adapter_info["points"]:
            for point in adapter_info["points"]:
                all_markdown_parts.append(f"- {point}\n")
        all_markdown_parts.append("\n---\n")  # Add a horizontal rule for separation

    all_markdown_parts.append(
        '<script>'
        '(function() {'
        '  function showTooltip(el) {'
        '    var tooltip = el.parentElement.querySelector(".scheme-raw-tooltip");'
        '    if (tooltip) {'
        '      tooltip.style.display = "block";'
        '      setTimeout(function() { tooltip.style.display = "none"; }, 1200);'
        '    }'
        '  }'
        '  document.querySelectorAll(".copy-scheme-raw").forEach(function(block) {'
        '    block.addEventListener("mouseenter", function() {'
        '      block.style.boxShadow = "0 0 0 2px #FF6F61";'
        '    });'
        '    block.addEventListener("mouseleave", function() {'
        '      block.style.boxShadow = "";'
        '    });'
        '    block.addEventListener("click", function(e) {'
        '      var scheme = block.getAttribute("data-scheme");'
        '      if (navigator.clipboard) {'
        '        navigator.clipboard.writeText(scheme).then(function() {'
        '          showTooltip(block);'
        '        });'
        '      } else {'
        '        var textarea = document.createElement("textarea");'
        '        textarea.value = scheme;'
        '        document.body.appendChild(textarea);'
        '        textarea.select();'
        '        document.execCommand("copy");'
        '        document.body.removeChild(textarea);'
        '        showTooltip(block);'
        '      }'
        '    });'
        '  });'
        '})();'
        '</script>'
    )

    try:
        adapters_md_path.parent.mkdir(parents=True, exist_ok=True)
        with open(adapters_md_path, "w", encoding="utf-8") as f:
            f.write("".join(all_markdown_parts))
        print(f"Successfully updated {adapters_md_path}")
    except Exception as e:
        print(f"Error writing to {adapters_md_path}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
