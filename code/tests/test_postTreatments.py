from unittest import TestCase
from parsing import get_data_from_csv
from postTreatments import PostTreatments
from fileCreation import load_targets_from_file
from config import *


class TestPostTreatment(TestCase):
    """
    Tests post treatments of targets except multiplex_machines and pseudo multiplex_machines\
    . tests with data1_10.csv

    **Huge disclaimer: tests will work with data1_10 AND BLAST = True(enable), overlap=False AND sort_by_tfgr = True**
    """
    def create_targets(self):
        """
        Creates and save targets for not re-execute the whole creation of primers pair
        execute only one time.
        """
        targets = get_data_from_csv(data1_10)
        for target in targets:
            target.compute_list_primers_pairs()
        return targets

    def test_check_overlap(self):
        """
        Tests if the overlapping is doing well
        """
        # execute one time the line commented below.
        self.create_targets()
        targets = load_targets_from_file(10)
        pt = PostTreatments(targets)
        res = pt.check_overlap()
        self.assertTrue(len(res) <= len(targets))
        for primer_pair in res:
            self.assertTrue(primer_pair.amplicon_range <= 90)

    def test_is_overlapping(self):
        """
        Tests if two primers_pair are overlapping each other.
        """
        targets = load_targets_from_file(10)
        pt = PostTreatments(targets)
        print(targets[0].all_primers_pairs[0])
        print(targets[1].all_primers_pairs[0])
        self.assertTrue(pt.is_overlapping(targets[0].all_primers_pairs[0], targets[1].all_primers_pairs[0]))

    def test_is_overlapping_merged(self):
        """
        Tests if two primers_pair are overlapping each other with one merged
        """
        targets = load_targets_from_file(10)
        pt = PostTreatments(targets)
        merged_target = pt.merge_primerspairs(targets[0].all_primers_pairs[0], targets[1].all_primers_pairs[0])
        self.assertTrue(pt.is_overlapping(merged_target, targets[2].all_primers_pairs[0]))

    def test_merge_primerspairs(self):
        """
        Tests if the merge of two primers pair are correct.
        """
        targets = load_targets_from_file(10)
        pt = PostTreatments(targets)
        merged_target = pt.merge_primerspairs(targets[0].all_primers_pairs[0], targets[1].all_primers_pairs[0])
        self.assertTrue(merged_target.TFGR == [117622179, 117622248])
        merged_target2 = pt.merge_primerspairs(merged_target, targets[2].all_primers_pairs[0])
        self.assertTrue(merged_target2.TFGR == [117622179, 117622248])


