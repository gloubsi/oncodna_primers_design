from unittest import TestCase
from parsing import get_data_from_csv
from config import *


class TestParser(TestCase):
    """
    Tests the csv parser
    """
    def test_get_data_from_csv_damaging(self):
        """
        Tests if a target with a damaging variant is correctly initialise
        """
        result = get_data_from_csv(data_21)
        self.assertTrue(len(result) == 1)
        target = result[0]
        self.assertTrue(target.no_chromosome == "chr10")
        self.assertTrue(target.mutation_pos == 89720807)
        self.assertTrue(target.var_freq == 0.2705)
        self.assertTrue(target.bio_imp == "Damaging")
        self.assertTrue(target.var_key == "CHR10_89720808_T_TA")

    def test_get_data_from_csv_polymorphism(self):
        """
        Tests if a target with a probably polymorphism variant is correctly initialise
        """
        result = get_data_from_csv(data_1)
        self.assertTrue(len(result) == 1)
        target = result[0]
        self.assertTrue(target.no_chromosome == "chr6")
        self.assertTrue(target.mutation_pos == 117622183)
        self.assertTrue(target.var_freq == 0.5105)
        self.assertTrue(target.bio_imp == "Probably Polymorphism")
        self.assertTrue(target.var_key == "CHR6_117622184_G_C")

    def test_get_data_from_csv_unknown(self):
        """
        Tests if a target with a unknown variant is correctly initialise
        """
        result = get_data_from_csv(data_16)
        self.assertTrue(len(result) == 1)
        target = result[0]
        self.assertTrue(target.no_chromosome == "chr11")
        self.assertTrue(target.mutation_pos == 314206)
        self.assertTrue(target.var_freq == 0.9159)
        self.assertTrue(target.bio_imp == "Unknown")
        self.assertTrue(target.var_key == "CHR11_314207_C_G")
