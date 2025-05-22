import logging
import re
import sys
from importlib import resources

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Python < 3.11

# This is the common module for cutseq
# It will contain shared classes and functions


def load_adapters() -> dict:
    """Loads adapter definitions from the 'adapters.toml' file."""
    try:
        # Python 3.9+ provides resources.files()
        adapter_info = tomllib.loads(
            resources.files("cutseq")
            .joinpath("adapters.toml")
            .read_text(encoding="utf-8")
        )
        # return only the name and scheme mapping
        return {k: v["scheme"] for k, v in adapter_info.items() if "scheme" in v}
    except Exception as e:
        logging.error(f"Error loading or parsing 'adapters.toml': {e}")
        # Depending on how critical this is, either raise e, sys.exit(), or return empty dict
        # For now, let's re-raise to make the error visible
        raise e


BUILDIN_ADAPTERS = load_adapters()


def reverse_complement(b):
    """
    Calculates the reverse complement of a DNA sequence.

    :param b: The input DNA sequence string.
    :type b: str
    :return: The reverse complemented DNA sequence string.
    :rtype: str
    """
    comp_map = dict(zip("ATGCatgc", "TACGtacg"))
    return "".join([comp_map.get(x, x) for x in b[::-1]])


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


class BarcodeSeq:
    """
    Represents a DNA sequence and its reverse complement.

    This class stores a forward DNA sequence, its reverse complement, and its length.
    It is used to conveniently access both forms of a sequence, typically for adapter
    or barcode sequences.

    :ivar fw: The forward DNA sequence.
    :vartype fw: str
    :ivar rc: The reverse complemented DNA sequence.
    :vartype rc: str
    :ivar len: The length of the forward sequence.
    :vartype len: int
    """

    def __init__(self, seq):
        """
        Initializes a BarcodeSeq object.

        :param seq: The forward DNA sequence.
        :type seq: str
        """
        self.fw = seq
        self.rc = reverse_complement(seq)
        self.len = len(seq)

    def __repr__(self):
        if self.len == 0:
            return ""
        return f"{self.fw} ({self.rc})"


class BarcodeConfig:
    """
    Parses and stores barcode and adapter scheme configurations.

    This class takes an adapter scheme string, parses it into its components
    (e.g., p5/p7 adapters, inline barcodes, UMIs, mask sequences, strand information),
    and stores these components as BarcodeSeq objects or strings.

    The expected adapter scheme format is a string like:
    `(p5_sequence)([inline5_sequence])(umi5_Ns)(mask5_Xs)(strand_indicator)(mask3_Xs)(umi3_Ns)([inline3_sequence])(p7_sequence)`
    Where:
        - `()` denote optional inline sequences.
        - `N` represents UMI bases.
        - `X` represents mask bases (to be trimmed).
        - `strand_indicator` is one of '>', '<', or '-' for forward, reverse, or unspecified strand.

    :ivar strand: Strand information ('+', '-', or None).
    :ivar p5: BarcodeSeq object for the P5 adapter.
    :ivar p7: BarcodeSeq object for the P7 adapter.
    :ivar inline5: BarcodeSeq object for the 5' inline barcode.
    :ivar inline3: BarcodeSeq object for the 3' inline barcode.
    :ivar umi5: BarcodeSeq object for the 5' UMI.
    :ivar umi3: BarcodeSeq object for the 3' UMI.
    :ivar mask5: BarcodeSeq object for the 5' mask sequence.
    :ivar mask3: BarcodeSeq object for the 3' mask sequence.
    """

    def __init__(self, adapter=None):
        """
        Initializes a BarcodeConfig object.

        :param adapter: The adapter scheme string to parse. If None, an empty
                        configuration is created.
        :type adapter: str, optional
        """
        self.strand = None
        self.p5 = BarcodeSeq("")
        self.p7 = BarcodeSeq("")
        self.inline5 = BarcodeSeq("")
        self.inline3 = BarcodeSeq("")
        self.umi5 = BarcodeSeq("")
        self.umi3 = BarcodeSeq("")
        self.mask5 = BarcodeSeq("")
        self.mask3 = BarcodeSeq("")
        if adapter is not None:
            self._parse_barcode(adapter)

    def _parse_barcode(self, b):
        """
        Parses the adapter scheme string.

        This method uses a regular expression to identify and extract the different
        components of the adapter scheme string. It populates the instance attributes
        based on the parsed values. If the scheme is invalid, it logs an error
        and exits.

        :param b: The adapter scheme string.
        :type b: str
        :raises SystemExit: If the barcode scheme string is invalid.
        """
        m = re.match(
            r"(?P<p5>[ATGCatgc]+)(\((?P<inline5>[ATGCatgc]+)\))?(?P<umi5>N*)(?P<mask5>X*)(?P<strand>-|>|<)(?P<mask3>X*)(?P<umi3>N*)(\((?P<inline3>[ATGCatgc]+)\))?(?P<p7>[ATGCatgc]+)",
            b,
        )
        if m is None:
            logging.error(f"barcode {b} is not valid")
            sys.exit(1)
        d = m.groupdict()
        if d["inline5"] is None:
            d["inline5"] = ""
        if d["inline3"] is None:
            d["inline3"] = ""
        self.strand = "+" if d["strand"] == ">" else "-" if d["strand"] == "<" else None
        self.p5 = BarcodeSeq(d["p5"])
        self.p7 = BarcodeSeq(d["p7"])
        self.inline5 = BarcodeSeq(d["inline5"])
        self.inline3 = BarcodeSeq(d["inline3"])
        self.umi5 = BarcodeSeq(d["umi5"])
        self.umi3 = BarcodeSeq(d["umi3"])
        self.mask5 = BarcodeSeq(d["mask5"])
        self.mask3 = BarcodeSeq(d["mask3"])

    def to_dict(self):
        """
        Converts the barcode configuration to a dictionary.

        :return: A dictionary representation of the barcode configuration,
                 containing the forward sequences of all components and the strand.
        :rtype: dict
        """
        return {
            "p5": self.p5.fw,
            "p7": self.p7.fw,
            "inline5": self.inline5.fw,
            "inline3": self.inline3.fw,
            "umi5": self.umi5.fw,
            "umi3": self.umi3.fw,
            "mask5": self.mask5.fw,
            "mask3": self.mask3.fw,
            "strand": self.strand,
        }


def print_builtin_adapters():
    """
    Print all built-in adapter names and their schemes in a well-organized, pretty table.
    """
    from textwrap import wrap

    print("\nBuilt-in adapter schemes:\n")
    # Find the max width for name and scheme for alignment
    max_name_len = max(len(name) for name in BUILDIN_ADAPTERS)
    max_scheme_len = max(len(scheme) for scheme in BUILDIN_ADAPTERS.values())
    # Header
    print(f"{'Name'.ljust(max_name_len)}   {'Scheme'}")
    print(f"{'-'*max_name_len}   {'-'*max(30, min(max_scheme_len, 100))}")
    # Print each adapter, wrapping long schemes
    for name, scheme in BUILDIN_ADAPTERS.items():
        wrapped_scheme = wrap(scheme, width=100)
        print(f"{name.ljust(max_name_len)}   {wrapped_scheme[0]}")
        for cont in wrapped_scheme[1:]:
            print(f"{' '*max_name_len}   {cont}")
    print("\nUse the adapter name with -A/--adapter-name, or the scheme string with -a/--adapter-scheme.\n")
