from unittest import TestCase
from parsing import get_data_from_csv
from blast import Blast
from config import *


class TestBlast(TestCase):
    """
    Tests the BLAST class
    """
    def test_get_unique_primers(self):
        """
        Tests if the function get correctly each unique primers and\
        each primer has not been forgotten.
        """
        target = get_data_from_csv(data_1)[0]
        target.compute_optimal_primers_pairs()
        blast = Blast(target.all_primers_pairs)
        u_l_primers, u_r_primers = blast.get_unique_primers()
        self.assertTrue(len(u_l_primers) == len(set(u_l_primers)))
        self.assertTrue(len(u_r_primers) == len(set(u_r_primers)))
        for primers_pair in target.all_primers_pairs:
            self.assertTrue(primers_pair.right_primer in u_r_primers)
            self.assertTrue(primers_pair.left_primer in u_l_primers)

    def test_run_blast(self):
        """
        Tests if blast hits are correct respect to our constraint.
        """
        target = get_data_from_csv(data_1)[0]
        target.compute_optimal_primers_pairs()
        blast = Blast(target.all_primers_pairs)
        u_l_primers, u_r_primers = blast.get_unique_primers()
        blast.run_blast(u_l_primers[0])
        with open(blast_output_file_path, "r") as o:
            for line in o:
                blast_result = line.split()
                self.assertTrue(float(blast_result[2]) == 100.)
                self.assertTrue(float(blast_result[-2]) < 1)
                self.assertTrue(int(blast_result[3])/float(blast_result[7]) >= 0.8)

        blast.run_blast(u_r_primers[0])
        with open(blast_output_file_path, "r") as o:
            for line in o:
                blast_result = line.split()
                self.assertTrue(float(blast_result[2]) == 100.)
                self.assertTrue(float(blast_result[-2]) < 1)
                self.assertTrue(int(blast_result[3]) / float(blast_result[7]) >= 0.8)

    def test_search_hybridization_site_rr_ll(self):
        """
        Tests if the function detects well rr or ll hybridisation
        Tests the function is_primers_compatible as well.
        """
        target = get_data_from_csv(data_16)[0]
        target.compute_optimal_primers_pairs()
        blast = Blast(target.all_primers_pairs)
        u_l_primers, u_r_primers = blast.get_unique_primers()
        blast.parse_blast_output(u_l_primers[0])
        self.assertFalse(blast.search_hybridization_site_rr_ll(u_l_primers[0]))
        u_l_primers[0].primer_hybridization_sites["chr20"].append([50283620, 50283640, 1])
        self.assertTrue(blast.search_hybridization_site_rr_ll(u_l_primers[0]))

        blast.parse_blast_output(u_r_primers[0])
        self.assertFalse(blast.search_hybridization_site_rr_ll(u_r_primers[0]))
        u_r_primers[0].primer_hybridization_sites["chr11"].append([60672820, 60672805, -1])
        self.assertTrue(blast.search_hybridization_site_rr_ll(u_r_primers[0]))

    def test_search_hybridization_sites_lr_rl(self):
        """
        Tests if the function detects well lr or rl hybridisation
        """
        target = get_data_from_csv(data_16)[0]
        target.compute_optimal_primers_pairs()
        blast = Blast(target.all_primers_pairs)
        u_l_primers, u_r_primers = blast.get_unique_primers()
        blast.parse_blast_output(u_l_primers[0])
        blast.parse_blast_output(u_r_primers[0])
        self.assertFalse(blast.search_hybridization_sites_lr_rl(u_l_primers[0].primer_hybridization_sites,
                                                                u_r_primers[0].primer_hybridization_sites))

        u_r_primers[0].primer_hybridization_sites["chr20"] = [[50283620, 50283640, 1]]
        self.assertTrue(blast.search_hybridization_sites_lr_rl(u_l_primers[0].primer_hybridization_sites,
                                                               u_r_primers[0].primer_hybridization_sites))

    def add_hybridization_sites(self):
        """
        Tests succinctly the main function of BLAST class
        """
        target = get_data_from_csv(data_16)[0]
        target.compute_optimal_primers_pairs()
        primers_pair = target.all_primers_pairs
        blast = Blast(target.all_primers_pairs)
        blast.add_hybridization_sites()
        print(len(primers_pair))

