import unittest
import sys
import os

# Adjust sys.path to find the cutseq module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from cutseq.common import BarcodeConfig, remove_fq_suffix, reverse_complement

class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement_valid(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("N"), "N")
        self.assertEqual(reverse_complement(""), "")
        self.assertEqual(reverse_complement("GATTACA"), "TGTAATC")

    def test_reverse_complement_mixed_case(self):
        self.assertEqual(reverse_complement("aTgC"), "gCaT")

class TestBarcodeConfig(unittest.TestCase):
    def test_parse_barcode_full_scheme(self):
        bc = BarcodeConfig("AGATCGGAAGAGCACACGTCTGACGC(TCGAC)NNNXXX>XXXNNN(GCTA)TGCATCTAG")
        self.assertEqual(bc.p5.fw, "AGATCGGAAGAGCACACGTCTGACGC")
        self.assertEqual(bc.inline5.fw, "TCGAC")
        self.assertEqual(bc.umi5.fw, "NNN")
        self.assertEqual(bc.mask5.fw, "XXX")
        self.assertEqual(bc.strand, "+")
        self.assertEqual(bc.mask3.fw, "XXX")
        self.assertEqual(bc.umi3.fw, "NNN")
        self.assertEqual(bc.inline3.fw, "GCTA")
        self.assertEqual(bc.p7.fw, "TGCATCTAG")

    def test_parse_barcode_minimal_scheme(self):
        bc = BarcodeConfig("AGATC>TGCAT")
        self.assertEqual(bc.p5.fw, "AGATC")
        self.assertEqual(bc.inline5.fw, "")
        self.assertEqual(bc.umi5.fw, "")
        self.assertEqual(bc.mask5.fw, "")
        self.assertEqual(bc.strand, "+")
        self.assertEqual(bc.mask3.fw, "")
        self.assertEqual(bc.umi3.fw, "")
        self.assertEqual(bc.inline3.fw, "")
        self.assertEqual(bc.p7.fw, "TGCAT")

    def test_parse_barcode_no_inline(self):
        bc = BarcodeConfig("AGATCNNNXXX>XXXNNNTGCAT")
        self.assertEqual(bc.p5.fw, "AGATC")
        self.assertEqual(bc.inline5.fw, "")
        self.assertEqual(bc.umi5.fw, "NNN")
        self.assertEqual(bc.mask5.fw, "XXX")
        self.assertEqual(bc.strand, "+")
        self.assertEqual(bc.mask3.fw, "XXX")
        self.assertEqual(bc.umi3.fw, "NNN")
        self.assertEqual(bc.inline3.fw, "")
        self.assertEqual(bc.p7.fw, "TGCAT")

    def test_parse_barcode_no_umi(self):
        bc = BarcodeConfig("AGATC(TCGA)XXX>XXX(GCTA)TGCAT")
        self.assertEqual(bc.p5.fw, "AGATC")
        self.assertEqual(bc.inline5.fw, "TCGA")
        self.assertEqual(bc.umi5.fw, "")
        self.assertEqual(bc.mask5.fw, "XXX")
        self.assertEqual(bc.strand, "+")
        self.assertEqual(bc.mask3.fw, "XXX")
        self.assertEqual(bc.umi3.fw, "")
        self.assertEqual(bc.inline3.fw, "GCTA")
        self.assertEqual(bc.p7.fw, "TGCAT")

    def test_parse_barcode_no_mask(self):
        bc = BarcodeConfig("AGATC(TCGA)NNN>NNN(GCTA)TGCAT")
        self.assertEqual(bc.p5.fw, "AGATC")
        self.assertEqual(bc.inline5.fw, "TCGA")
        self.assertEqual(bc.umi5.fw, "NNN")
        self.assertEqual(bc.mask5.fw, "")
        self.assertEqual(bc.strand, "+")
        self.assertEqual(bc.mask3.fw, "")
        self.assertEqual(bc.umi3.fw, "NNN")
        self.assertEqual(bc.inline3.fw, "GCTA")
        self.assertEqual(bc.p7.fw, "TGCAT")

    def test_parse_barcode_strand_minus(self):
        bc = BarcodeConfig("AGATC<TGCAT")
        self.assertEqual(bc.p5.fw, "AGATC")
        self.assertEqual(bc.strand, "-")
        self.assertEqual(bc.p7.fw, "TGCAT")

    def test_parse_barcode_strand_none(self):
        bc = BarcodeConfig("AGATC-TGCAT")
        self.assertEqual(bc.p5.fw, "AGATC")
        self.assertEqual(bc.strand, None)
        self.assertEqual(bc.p7.fw, "TGCAT")

    def test_parse_barcode_only_p5_p7_strand(self):
        bc = BarcodeConfig("GATTACA>TGTAATC")
        self.assertEqual(bc.p5.fw, "GATTACA")
        self.assertEqual(bc.inline5.fw, "")
        self.assertEqual(bc.umi5.fw, "")
        self.assertEqual(bc.mask5.fw, "")
        self.assertEqual(bc.strand, "+")
        self.assertEqual(bc.mask3.fw, "")
        self.assertEqual(bc.umi3.fw, "")
        self.assertEqual(bc.inline3.fw, "")
        self.assertEqual(bc.p7.fw, "TGTAATC")
        self.assertEqual(bc.p5.rc, "TGTAATC") # testing reverse_complement via BarcodeSeq
        self.assertEqual(bc.p7.rc, "GATTACA") # testing reverse_complement via BarcodeSeq

    # TODO: Add tests for invalid schemes if sys.exit can be mocked or refactored.
    # For now, we rely on the BarcodeConfig's own error logging for invalid inputs.
    # Example of what an invalid test might look like (would require mocking sys.exit):
    # @unittest.mock.patch('sys.exit')
    # @unittest.mock.patch('logging.error')
    # def test_parse_barcode_invalid_scheme_char(self, mock_log_error, mock_sys_exit):
    #     BarcodeConfig("AGATC*TGCAT") # Invalid character *
    #     mock_log_error.assert_called_once()
    #     mock_sys_exit.assert_called_once_with(1)

class TestRemoveFqSuffix(unittest.TestCase):
    def test_remove_common_suffixes(self):
        self.assertEqual(remove_fq_suffix("myfile_R1.fastq.gz"), "myfile")
        self.assertEqual(remove_fq_suffix("myfile_R2.fq.gz"), "myfile")
        self.assertEqual(remove_fq_suffix("myfile_R1_001.fq"), "myfile")
        self.assertEqual(remove_fq_suffix("myfile.fastq"), "myfile")
        self.assertEqual(remove_fq_suffix("myfile_R1_001.fastq.gz"), "myfile") # Longest suffix first
        self.assertEqual(remove_fq_suffix("myfile_R1.fq"), "myfile")
        self.assertEqual(remove_fq_suffix("myfile_R2_001.fastq"), "myfile")

    def test_remove_no_suffix(self):
        self.assertEqual(remove_fq_suffix("myfile"), "myfile")
        self.assertEqual(remove_fq_suffix("myfile.txt"), "myfile.txt")
        self.assertEqual(remove_fq_suffix("myfile_R3.fastq.gz"), "myfile_R3.fastq.gz") # Not a standard R1/R2

    def test_remove_suffix_with_dots_in_name(self):
        self.assertEqual(remove_fq_suffix("my.file.name_R1.fastq.gz"), "my.file.name")

    def test_remove_empty_string(self):
        self.assertEqual(remove_fq_suffix(""), "")

    def test_ensure_longest_suffix_is_removed(self):
        # _R1_001.fastq.gz vs .fastq.gz
        self.assertEqual(remove_fq_suffix("test_R1_001.fastq.gz"), "test")
        # _R1.fastq.gz vs .fastq.gz
        self.assertEqual(remove_fq_suffix("test_R1.fastq.gz"), "test")
        # _R1_001.fq vs .fq
        self.assertEqual(remove_fq_suffix("test_R1_001.fq"), "test")
        # .fq.gz vs .gz (though .gz alone is not a suffix we remove, this ensures .fq.gz is caught)
        self.assertEqual(remove_fq_suffix("test.fq.gz"), "test")


if __name__ == '__main__':
    unittest.main()
