import logging
from importlib import resources

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Python < 3.11

# This is the common module for cutseq
# It will contain shared classes and functions


def _load_adapter_meta() -> dict:
    """Load the full adapters.toml (all sections + fields)."""
    try:
        return tomllib.loads(
            resources.files("cutseq")
            .joinpath("adapters.toml")
            .read_text(encoding="utf-8")
        )
    except Exception as e:
        logging.error(f"Error loading or parsing 'adapters.toml': {e}")
        raise e


def load_adapters() -> dict:
    """Loads adapter definitions from the 'adapters.toml' file.

    Each ``[NAME]`` section with a ``scheme`` becomes a built-in name.
    A section with only ``alias_of = "TARGET"`` (and no scheme) resolves to
    the target's scheme, so legacy names keep working after a rename.
    """
    adapter_info = _load_adapter_meta()
    out = {}
    for k, v in adapter_info.items():
        if not isinstance(v, dict):
            continue
        if "scheme" in v:
            out[k] = v["scheme"]
        elif v.get("alias_of"):
            out[k] = adapter_info.get(v["alias_of"], {}).get("scheme")
    return out


BUILDIN_ADAPTERS = load_adapters()


def remove_fq_suffix(f):
    """
    Removes common FASTQ file extensions from a filename.

    This function attempts to remove suffixes like '_R1.fastq.gz', '.fq', etc.,
    to derive a base filename. It prioritizes longer, more specific suffixes.

    :param f: The input filename string.
    :type f: str
    :return: The filename with standard FASTQ suffixes removed.
    :rtype: str
    :Example:
        >>> remove_fq_suffix("my_sample_R1.fastq.gz")
        'my_sample'
        >>> remove_fq_suffix("another_file.fq")
        'another_file'
        >>> remove_fq_suffix("no_suffix_here")
        'no_suffix_here'
    """
    suffixes = [
        f"{base}.{ext}"
        for ext in ["fastq.gz", "fq.gz", "fastq", "fq"]
        for base in ["_R1_001", "_R2_001", "_R1", "_R2", ""]
    ]

    for suffix in suffixes:
        if f.endswith(suffix):
            return f.removesuffix(suffix)
    return f



def print_builtin_adapters():
    """
    Print all built-in adapter names and their schemes in a well-organized, pretty table.
    """
    from textwrap import wrap

    meta = _load_adapter_meta()
    labels = {}
    for name in BUILDIN_ADAPTERS:
        entry = meta.get(name, {})
        if entry.get("alias_of"):
            labels[name] = f"{name}  (alias -> {entry['alias_of']})"
        else:
            labels[name] = name

    print("\nBuilt-in adapter schemes:\n")
    # Find the max width for name and scheme for alignment
    max_name_len = max(len(labels[n]) for n in BUILDIN_ADAPTERS)
    max_scheme_len = max(len(scheme) for scheme in BUILDIN_ADAPTERS.values())
    # Header
    print(f"{'Name'.ljust(max_name_len)}   {'Scheme'}")
    print(f"{'-'*max_name_len}   {'-'*max(30, min(max_scheme_len, 100))}")
    # Print each adapter, wrapping long schemes
    for name, scheme in BUILDIN_ADAPTERS.items():
        wrapped_scheme = wrap(scheme, width=100)
        print(f"{labels[name].ljust(max_name_len)}   {wrapped_scheme[0]}")
        for cont in wrapped_scheme[1:]:
            print(f"{' '*max_name_len}   {cont}")
    print("\nUse the adapter name or a custom grammar scheme with -A/--adapter-scheme (-a is an alias).\n")
