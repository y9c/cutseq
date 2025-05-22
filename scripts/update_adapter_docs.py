#!/usr/bin/env python3
import sys
from pathlib import Path

# Add project root to sys.path to allow importing cutseq.common
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

try:
    from cutseq.common import BarcodeConfig
except ImportError as e:
    print(f"Error importing BarcodeConfig from cutseq.common: {e}")
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
    "p5": "#A8E6CF",
    "p7": "#D1E8D1",
    "umi": "#B2EBF2",
    "inline": "#FFD700",
    "mask": "#DCDCDC",
    "strand_bg": "#FF6F61",  # A distinct background for the strand indicator
    "default_seq": "#E0E0E0",
}


def get_part_type(part_name):
    """Determines the 'type' of adapter part for color selection."""
    if "p5" in part_name:
        return "p5"
    if "p7" in part_name:
        return "p7"
    if "umi" in part_name:
        return "umi"
    if "inline" in part_name:
        return "inline"
    if "mask" in part_name:
        return "mask"
    return "default_seq"  # Fallback


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
    all_markdown_parts.append("---\ntitle: Adapter Schemes\nnav_order: 2\n---\n\n")
    all_markdown_parts.append("# Adapter Schemes\n\n")
    all_markdown_parts.append(
        "CutSeq supports a variety of built-in adapter schemes for common NGS library types. You can list all available schemes in your terminal with:\n\n"
    )
    all_markdown_parts.append(
        "```bash\ncutseq --list-adapters\n```\n\n"
    )
    all_markdown_parts.append(
        "Use the adapter name with `-A/--adapter-name`, or specify a custom scheme string with `-a/--adapter-scheme`.\n\n"
    )
    all_markdown_parts.append(
        "## Example: Built-in Schemes\n\n"
    )
    all_markdown_parts.append(
        "- **SMALLRNA**: Small RNA libraries, double ligation, forward orientation\n"
        "- **INLINE**: Custom barcoded libraries, dual UMI, inline barcode\n"
        "- **TAKARAV3**: SMARTer Stranded Total RNA-Seq Kit v3\n"
        "- **STRANDED**: Stranded RNA libraries\n\n"
    )
    all_markdown_parts.append(
        "See below for a comprehensive guide to each supported adapter pattern, including copyable scheme blocks and usage notes.\n\n---\n\n"
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

        # HTML Visualization
        bc = BarcodeConfig(adapter_info["scheme"])
        scheme_raw = adapter_info["scheme"]

        html_parts = [
            '<div class="adapter-scheme" style="margin-bottom: 15px; position: relative;">',
            f'<div class="copy-scheme-raw" style="display: flex; flex-wrap: nowrap; align-items: center; font-family: monospace; font-size: 14px; border: 1px solid #ccc; padding: 5px; border-radius: 5px; overflow-x: auto; cursor: pointer; background: #f8f8f8; transition: box-shadow 0.2s;" title="Click to copy scheme: {scheme_raw}" data-scheme="{scheme_raw}">'  # clickable block
        ]

        adapter_components_ordered = [
            "p5",
            "inline5",
            "umi5",
            "mask5",
            "strand",
            "mask3",
            "umi3",
            "inline3",
            "p7",
        ]

        for part_name in adapter_components_ordered:
            seq_val = ""
            is_strand_part = False

            if part_name == "strand":
                strand_char_map = {"+": ">", "-": "<", None: "-"}
                seq_val = strand_char_map.get(getattr(bc, 'strand', None), "-")
                is_strand_part = True
            elif hasattr(bc, part_name):
                barcode_seq_obj = getattr(bc, part_name)
                if hasattr(barcode_seq_obj, "fw"):
                    seq_val = barcode_seq_obj.fw

            if seq_val:
                if is_strand_part:
                    html_parts.append(
                        f'<div style="position: relative; width: 30px; height: 30px; margin: 0 2px; text-align: center; line-height: 30px;"><div style="background-color: {COLOR_PALETTE["strand_bg"]}; width: 100%; height: 100%; position: absolute; top: 0; left: 0; clip-path: polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%);"></div><span style="position: relative; z-index: 1; color: white; font-weight: bold;">{seq_val}</span></div>'
                    )
                else:
                    part_type = get_part_type(part_name)
                    color = COLOR_PALETTE.get(part_type, COLOR_PALETTE["default_seq"])
                    html_parts.append(
                        f'<span style="background-color: {color}; padding: 5px 8px; margin: 0 2px; border-radius: 3px; white-space: nowrap;">{seq_val}</span>'
                    )

        html_parts.append("</div>")  # Close flex container
        # Add a hidden raw scheme and a tooltip for copy feedback
        html_parts.append(
            '<div class="scheme-raw-tooltip" style="display:none; position:absolute; top:-30px; left:0; background:#222; color:#fff; padding:3px 8px; border-radius:4px; font-size:12px; z-index:10;">Copied!</div>'
        )
        html_parts.append("</div>\n")  # Close adapter-scheme div

        all_markdown_parts.extend(html_parts)

        # Bullet Points
        if adapter_info["points"]:
            for point in adapter_info["points"]:
                all_markdown_parts.append(f"- {point}\n")
        all_markdown_parts.append("\n---\n")  # Add a horizontal rule for separation

    # --- End Generate Markdown Content ---
    # Remove unused markdown_content assignment
    # Add JS/CSS for copy-to-clipboard functionality at the end of the markdown
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
