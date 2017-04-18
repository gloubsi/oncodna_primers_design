from unittest import TestCase
from parsing import get_data_from_csv
from config import *


class TestTarget(TestCase):
    """
    Tests the class Target especially creation of primers by primer3 algorithm
    the compute_list_primers_pairs function has not been tested because it's a wrapper of good tested functions.
    """
    def test_get_target_sequence(self):
        """
        Tests if the correct target sequence is grabbed
        """
        target = get_data_from_csv(data_1)[0]
        target.get_target_sequence(100)
        self.assertTrue(len(target.sequence) == 201)

    def test_design_primers(self):
        """
        Tests if designed primers pair are compatible with constraints
        """
        target = get_data_from_csv(data_1)[0]
        target.get_target_sequence(100)
        primers_pair = target.design_primers(target.range, target.sequence)
        print(primers_pair)
        self.primers_caract(primers_pair, target.range)

    def test_compute_optimal_raw_primers_pairs(self):
        """
        Tests for all primers pair generated that the primers pairs respect constraints
        """
        def test_compute_optimal_raw_primers_pairs_helper(file_data_path):
            target = get_data_from_csv(file_data_path)[0]
            result = target.compute_optimal_raw_primers_pairs()
            raw_primers_pairs = result[0]
            ranges_left = result[1]
            for i in range(len(raw_primers_pairs)):
                variant_position = ranges_left[i] + 1
                self.primers_caract(raw_primers_pairs[i], variant_position)

        test_compute_optimal_raw_primers_pairs_helper(data_1)
        test_compute_optimal_raw_primers_pairs_helper(data_16)
        test_compute_optimal_raw_primers_pairs_helper(data_21)

    def primers_caract(self, primers_pair, variant_position):
        """
        Tests main characteristics of primers pair.
        :param primers_pair: a primer pair (not the primersPair instance) as primer3 output.
        :param variant_position: the relative position of the variant
        """
        # primers pair not null
        self.assertTrue("PRIMER_LEFT_0_SEQUENCE" in primers_pair)
        self.assertTrue("PRIMER_RIGHT_0_SEQUENCE" in primers_pair)

        # check if primer around variant
        self.assertTrue(primers_pair["PRIMER_LEFT_0"][0] + primers_pair["PRIMER_LEFT_0"][1] < variant_position <
                        primers_pair["PRIMER_RIGHT_0"][0] - primers_pair["PRIMER_RIGHT_0"][1])

        # GC content
        self.assertTrue(30 <= primers_pair["PRIMER_LEFT_0_GC_PERCENT"] <= 70)
        self.assertTrue(30 <= primers_pair["PRIMER_RIGHT_0_GC_PERCENT"] <= 70)

        # Melting temperature (TM)
        self.assertTrue(57 <= primers_pair["PRIMER_LEFT_0_TM"] <= 63)
        self.assertTrue(57 <= primers_pair["PRIMER_RIGHT_0_TM"] <= 63)

    def test_sort_by_tfgr(self):
        """
        Tests if primers pairs instance have been correctly sorted by TFGR.
        """
        target = get_data_from_csv(data_1)[0]
        target.compute_optimal_primers_pairs()
        target.sort_by_tfgr()
        amplicon_range = target.all_primers_pairs[0].amplicon_range
        for primers_pair in target.all_primers_pairs[1:]:
            self.assertTrue(primers_pair.amplicon_range >= amplicon_range)
            amplicon_range = primers_pair.amplicon_range

    def test_sort_by_penalty(self):
        """
        Tests if primers pairs instance have been correctly sorted by penalty.
        """
        target = get_data_from_csv(data_1)[0]
        target.compute_optimal_primers_pairs()
        target.sort_by_penalty()
        global_penalty = target.all_primers_pairs[0].penalty
        for primers_pair in target.all_primers_pairs[1:]:
            self.assertTrue(primers_pair.penalty >= global_penalty)
            global_penalty = primers_pair.penalty

    def test_remove_primers_pair_with_hybridisation(self):
        """
        Tests if the function remove all primers_pair with any hybridisation sites.
        """
        target = get_data_from_csv(data_1)[0]
        target.compute_optimal_primers_pairs()
        target.remove_primers_pair_with_hybridisation()
        for primers_pair in target.all_primers_pairs:
            self.assertTrue(not primers_pair.left_primer.primer_hybridization_sites)
            self.assertTrue(not primers_pair.right_primer.primer_hybridization_sites)

