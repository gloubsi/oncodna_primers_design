from unittest import TestCase
from parsing import get_data_from_csv
from primer import Primer
from primersPair import PrimersPair
from config import *


class TestPrimersPair(TestCase):
    """
    Tests if the primersPair instance is correctly created
    """
    def test_eq(self):
        target = get_data_from_csv(data_1)[0]
        primers_pair1 = self.create_primers_pair(target)
        primers_pair2 = self.create_primers_pair(target)
        self.assertTrue(primers_pair1 == primers_pair2)

        target2 = get_data_from_csv(data_16)[0]
        primers_pair3 = self.create_primers_pair(target2)
        self.assertTrue(primers_pair3 != primers_pair2)
        l = [primers_pair1]
        self.assertTrue(primers_pair1 in l)
        self.assertTrue(primers_pair2 in l)
        self.assertTrue(primers_pair3 not in l)

    def create_primers_pair(self, target):
        """
        Creates a primers_pair from a target instance

        :param target: target instance

        :return: a primers_pair instance (in target sequence)
        """
        target.get_target_sequence(target.range)
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

        return PrimersPair(left_primer, right_primer)

