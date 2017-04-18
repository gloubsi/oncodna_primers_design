#!/usr/local/bin/python

from fileCreation import create_fasta_file_primer
import subprocess
import shlex
from config import *


class Blast:
    """
    Class that computes for each primers pairs of a target possible hybridisation sites of a primer

    :param list_primers_pairs: a list containing primersPair instance.
    """
    def __init__(self, list_primers_pairs):
        self.all_primers_pairs = list_primers_pairs

    def add_hybridization_sites(self):
        """
        Main function that adding hybridisation site to primer object with python dic.

        1. Gets all unique left and right primers of primers pair list of a target.

        2. Execute blast (via get_blast_output()) for each unique left primers and right primers\
           remove any primers that hybridise other regions elsewhere in the genome.

        3. Parse again the list of all primers pairs of the target and assign the hybridisation site of\
           left and right primer to the sames ones (a lot of doublon in primers pairs list).\
           So a lot of blast queries are avoid using this method.

        4. search possible hybridisation site elsewhere in the genome for each primers pairs.\
           if one, remove the primers pair.

        NO check for hybridisation with others primers_pair, no multiplex_machines

        Blast give possible always in 5' -> 3' (direction of the elongation)
        There are 4 cases of possible hybridisation:

        * L-R: first value < second value for the left; first value > second value for the right

        * R-L: first value < second value for the right; first value > second value for the left

        * L-L: first value < second value for left; first value > second value for the left

        * R-R: first value < second value for right; first value > second value for the right

        :return: nothing but update all_primers_pairs that contained all primers pairs that not have hybridisation\
                 elsewhere in the genome.
        """
        # get unique primers
        left_primers, right_primers = self.get_unique_primers()
        # Assign unique left primer and right primer to primers pair of the target
        for primers_pair in self.all_primers_pairs:
            pair_removed = False
            for lp in left_primers:
                if primers_pair.left_primer == lp:
                    self.parse_blast_output(lp)
                    if self.search_hybridization_site_rr_ll(lp):
                        self.all_primers_pairs.remove(primers_pair)
                        pair_removed = True
                    else:
                        primers_pair.left_primer.primer_hybridization_sites = lp.primer_hybridization_sites
                    break

            for rp in right_primers:
                if not pair_removed:
                    if primers_pair.right_primer == rp:
                        self.parse_blast_output(rp)
                        if self.search_hybridization_site_rr_ll(rp):
                            self.all_primers_pairs.remove(primers_pair)
                            pair_removed = True
                        else:
                            primers_pair.right_primer.primer_hybridization_sites = rp.primer_hybridization_sites
                        break
                else:
                    break

            if not pair_removed:
                if self.search_hybridization_sites_lr_rl(primers_pair.left_primer.primer_hybridization_sites,
                                                          primers_pair.right_primer.primer_hybridization_sites):

                    self.all_primers_pairs.remove(primers_pair)

    def elongation_direction(self, pos1, pos2):
        """
        gets the elongation direction of a primer defined by pos1 and pos2 (extremity positions)

        :param pos1: the first position given by BLAST
        :param pos2: the second position given by BLAST

        :return: * 1 if the elongation is performed left to right
                 * 0 if the elongation is performed right to left
        """
        if pos1 < pos2:
            # elongation left to right
            return 1
        elif pos2 < pos1:
            # elongation right to left, so pos1 > pos2 !
            return -1

    def get_unique_primers(self):
        """
        output: list containing all left primers_pair and list containing all right primers_pair. No doublon accepted.
        [[primer_sequence, primer_position],...]
        Be careful, the primer_position is the last position of the left primer or the first position of the right one.
        Parse each pair of accepted primers_pair and add the left primer and the right primer in different list

        :return: a tuple containing a list of all unique left primers and all unique right primers.
        """
        left_primers = []
        right_primers = []
        for primers_pair in self.all_primers_pairs:
            if primers_pair.left_primer not in left_primers:
                left_primers.append(primers_pair.left_primer)
            if primers_pair.right_primer not in right_primers:
                right_primers.append(primers_pair.right_primer)
        return left_primers, right_primers

    def run_blast(self, primer):
        """
        Blast the primer against human genome (fasta file hg19 as reference). Creates the output file 'outputBlast.txt\
        fields of the output format: query acc., subject acc., identity pourcentage, alignment length, mismatches,\
        gap opens, q. start, q. end, s. start, s. end, e-value, bit score\
        blast-short type command.

        :param primer: primer to BLAST (Primer instance)

        :return: nothing but ** create the file outputBlast.txt ** containing all results from BLAST.
        """
        create_fasta_file_primer(primer.sequence)
        f = open(blast_output_file_path, 'w')
        cmd = "{executable} -task blastn-short -db {databaseFile} " \
              "-query {primersFile} -outfmt 6 " \
              "-evalue {e_value} -perc_identity {perc_i} -qcov_hsp_perc {gcov}" \
            .format(executable=blast_executable, databaseFile=genome_file_path, primersFile=primer_to_blast_file_path,
                    e_value=e_value_limit, perc_i=perc_identity, gcov=p_align)
        subprocess.call(shlex.split(cmd), shell=False, stdout=f)
        f.close()

    def parse_blast_output(self, primer):
        """
        run and parse the result and create a dic with chromosome number as key and list of 'content'
        as value content = [BLAST_pos1, BLAST_pos2, elongation_direction]

        :param primer: the primer that we want to blast his sequence.

        :return: nothing but have a side effect on primer hybridisation attribute.
        """
        self.run_blast(primer)
        primer_hybridization_sites = {}
        with open(blast_output_file_path, "r") as o:
            for line in o:
                blast_result = line.split()
                if int(blast_result[8]) != primer.TFGP + 1 and int(blast_result[9]) != primer.TFGP + 1:
                    content = [int(blast_result[8]), int(blast_result[9]),
                               self.elongation_direction(blast_result[8], blast_result[9])]
                    if blast_result[1] in primer_hybridization_sites:
                        primer.primer_hybridization_sites[blast_result[1]].append(content)
                    else:
                        primer.primer_hybridization_sites[blast_result[1]] = [content]

    def search_hybridization_site_rr_ll(self, primer):
        """
        Search RR or LL hybridisation for one primer.

        :param primer: Primer instance

        :return: * True if hybridisation

                 * False otherwise
        """
        for k in primer.primer_hybridization_sites.keys():
            l = primer.primer_hybridization_sites.get(k)
            if len(l) > 1:
                for i in range(len(l) - 1):
                    for j in range(i, len(l)):
                        if i != j:
                            if self.is_primers_compatible(l[i], l[j]):
                                return True
        return False

    def search_hybridization_sites_lr_rl(self, d_left, d_right):
        """
        :param d_left: the dic of hybridisation sites of left primer.
        :param d_right: the dic of hybridisation sites of right primer.

        :return: * False if no hybridisation elsewhere in the genome (check OK)

                 * True otherwise (hybridisation check no ok)
        """
        for k1 in d_right.keys():
            for k2 in d_left.keys():
                if k1 == k2:
                    for blastResult1 in d_right.get(k1):
                        for blastResult2 in d_left.get(k2):
                            if self.is_primers_compatible(blastResult1, blastResult2):
                                return True
                    break
        return False

    def is_primers_compatible(self, blast_result1, blast_result2):
        """
        First checks the direction of elongation. If one primer is in the 'plus' strand and the other in the 'minus' strand
        Then check if the distance is smaller than 1000 pb.

        :param blast_result1: 'content' list for first primer
        :param blast_result2: 'content' list for the second primer

        :return: * True if hybridization

                 * False otherwise
        """
        if blast_result1[2] == 1 and blast_result2[2] == -1:
            if 0 < blast_result2[1] - blast_result1[1] < 1000:
                return True
        elif blast_result1[2] == -1 and blast_result2[2] == 1:
            if 0 < blast_result1[1] - blast_result2[1] < 1000:
                return True
        return False
