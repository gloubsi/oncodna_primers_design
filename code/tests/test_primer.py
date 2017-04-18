from unittest import TestCase
from parsing import get_data_from_csv
from primer import Primer
from config import *


class TestPrimer(TestCase):
    """
    Tests the Primer class
    """
    def test_eq(self):
        """
        Tests equality and inequality between two primers
        """
        target = get_data_from_csv(data_1)[0]
        target.get_target_sequence(target.range)
        primers_pair = target.design_primers(target.range, target.sequence)

        left_primer = Primer(primers_pair["PRIMER_LEFT_0_SEQUENCE"],
                             "left",
                             primers_pair["PRIMER_LEFT_0"][0],
                             100,
                             target)
        right_primer = Primer(primers_pair["PRIMER_RIGHT_0_SEQUENCE"],
                              "right",
                              primers_pair["PRIMER_RIGHT_0"][0],
                              100,
                              target)

        self.assertTrue(right_primer == right_primer)
        self.assertTrue(left_primer == left_primer)
        self.assertTrue(right_primer != left_primer)

    def test_TFGR(self):
        """
        Tests if the TFGP is correctly computed
        """
        self.compute_thermofisher_good_pos(data_1)
        self.compute_thermofisher_good_pos(data_16)
        self.compute_thermofisher_good_pos(data_1, 70, 80)
        self.compute_thermofisher_good_pos(data_16, 40, 40)

    def compute_thermofisher_good_pos(self, file_data_path, range_left=100, range_right=100):
        target = get_data_from_csv(file_data_path)[0]
        target.get_target_sequence(target.range)
        sequence = ""
        if range_right == 100:
            sequence = target.sequence[(target.range - range_left):]
        else:
            sequence = target.sequence[(target.range - range_left):-(target.range - range_right)]
        primers_pair = target.design_primers(range_left, sequence)
        left_primer = Primer(primers_pair["PRIMER_LEFT_0_SEQUENCE"],
                             "left",
                             primers_pair["PRIMER_LEFT_0"][0],
                             range_left,
                             target)
        right_primer = Primer(primers_pair["PRIMER_RIGHT_0_SEQUENCE"],
                              "right",
                              primers_pair["PRIMER_RIGHT_0"][0],
                              range_left,
                              target)

        self.assertTrue(left_primer.TFGP < target.mutation_pos < right_primer.TFGP)
        self.assertTrue(target.mutation_pos - left_primer.TFGP == range_left + 1 - (primers_pair["PRIMER_LEFT_0"][0] + primers_pair["PRIMER_LEFT_0"][1]))
        self.assertTrue(right_primer.TFGP - target.mutation_pos == (primers_pair["PRIMER_RIGHT_0"][0] - primers_pair["PRIMER_RIGHT_0"][1]) - range_left + 1)

    def test_check_snp_in_primer(self):
        """
        Tests snp and no snp in primer sequence
        """
        target = get_data_from_csv(data_1)[0]
        target.get_target_sequence(target.range)
        snp = target.get_snp()
        primers_pair = target.design_primers(target.range, target.sequence)
        left_primer = Primer(primers_pair["PRIMER_LEFT_0_SEQUENCE"],
                             "left",
                             primers_pair["PRIMER_LEFT_0"][0],
                             target.range,
                             target)
        right_primer = Primer(primers_pair["PRIMER_RIGHT_0_SEQUENCE"],
                              "right",
                              primers_pair["PRIMER_RIGHT_0"][0],
                              target.range,
                              target)

        self.assertTrue(left_primer.check_snp_in_primer(snp) == 0)
        self.assertTrue(right_primer.check_snp_in_primer(snp) == 0)

    def test_dinucleotides_repeats(self):
        """
        Tests if the function detects di-nucleotides repeats
        """
        target = get_data_from_csv(data_1)[0]
        target.get_target_sequence(target.range)
        primers_pair = target.design_primers(target.range, target.sequence)
        left_primer = Primer(primers_pair["PRIMER_LEFT_0_SEQUENCE"],
                             "left",
                             primers_pair["PRIMER_LEFT_0"][0],
                             target.range,
                             target)
        self.assertFalse(left_primer.dinucleotides_repeats())
        left_primer.sequence = "GTTCAGTT" + "GCGCGCGCGC" + "CACAGTGCAGCG"
        self.assertTrue(left_primer.dinucleotides_repeats())
