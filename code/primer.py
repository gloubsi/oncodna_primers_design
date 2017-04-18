#!/usr/local/bin/python
from config import *


class Primer:
    """
    Primer is an object representing a primer either left or right.

    Two primers are equal if their sequence are the same and their TFGP are equal.

    :param target: the target instance where the primer come from.
    :param sequence: sequence of the primer
    :param left_or_right: left if it's a left primer, right otherwise
    :param relative_pos: relative position (get from compute_optimal_primers_pairs() from target.py)
    :param range_left: get from compute_optimal_primers_pairs() from target.py
    :param penalty: describes how much differ the characteristic of this primer upon optimal characteristics.

    :ivar TFGP: the last position of the primer if left\
                or the first position of the primer if right.

    :ivar primer_hybridization_sites: a dic with chromosome number as key and a list of hybridisation sites as value.
    """
    def __init__(self, sequence, left_or_right, relative_pos, range_left, target, penalty=0):
        self.target = target
        self.sequence = sequence
        self.left_or_right = left_or_right
        self.relative_pos = relative_pos
        self.range_left = range_left
        self.penalty = penalty
        self.TFGP = self.compute_thermofisher_good_pos()
        self.primer_hybridization_sites = {}
    
    def __str__(self):
        return "sequence: " + self.sequence +\
               "\nleft or right? " + self.left_or_right + \
               "\nThermoFisher Good Position: " + str(self.TFGP) + \
               "\nhybridization sites: " + str(self.primer_hybridization_sites)
    
    def __eq__(self, other):
        """
        Two primers are equal if their sequence are the same and their TFGP are equal.

        :param other: an another primer instance.

        :return: True if the two primers instance are equal
                 False otherwise
        """
        return (self.sequence == other.sequence) and (self.TFGP == other.TFGP)

    def __hash__(self):
        return hash(str(self))
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def reverse_complement(self):
        """
        :return: reverse complement of the primer sequence
        """
        reverse_seq = self.sequence[::-1]
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
        res = ''
        for base in reverse_seq:
            res += complement[base]
        return res
    
    def compute_thermofisher_good_pos(self):
        """
        :return: the last position of the left primer or the first position of the right one
        """
        if self.left_or_right == "left":
            first_pos = self.relative_pos
            last_pos = int(first_pos) + len(self.sequence)
            scaled_pos = self.range_left + 1 - last_pos
            pos_primer_scaled = int(self.target.mutation_pos) - int(scaled_pos)
            return pos_primer_scaled
        else:
            last_pos = self.relative_pos
            first_pos = int(last_pos) - len(self.sequence)
            scaled_pos = int(first_pos) - self.range_left + 1
            pos_primer_scaled = int(self.target.mutation_pos) + int(scaled_pos)
            return pos_primer_scaled
            
    def check_snp_in_primer(self, snp_in_target):
        """
        Checks if a snp is in the primer

        :param snp_in_target: list containing all snp of the target\
                              get by get_snp() from target.py

        :return: * -1 if a snp is in the primer sequence

                 * 0 otherwise
        """
        if self.left_or_right == "left":
            primer_interval = [self.TFGP - len(self.sequence), self.TFGP]
        else:
            primer_interval = [self.TFGP, self.TFGP + len(self.sequence)]
        for element in snp_in_target:
            if primer_interval[0] < element[1] < primer_interval[1] and element[2] > gmaf_limit:
                return -1
        return 0
    
    def dinucleotides_repeats(self):
        """
        Checks if there is 5 or more consecutive dinucleotides repeats in the primer sequence.
        No more than 4 consecutive repeats.

        :return: * True if the sequence containing such repeats

                 * False otherwise
        """
        possible_repeats = ["ATATATATAT", "ACACACACAC", "AGAGAGAGAG",
                            "TATATATATA", "TCTCTCTCTC", "TGTGTGTGTG",
                            "GAGAGAGAGA", "GCGCGCGCGC", "GTGTGTGTGT",
                            "CGCGCGCGCG", "CACACACACA", "CTCTCTCTCT"]
        for repeat in possible_repeats:
            if repeat in self.sequence:
                return True
        return False
