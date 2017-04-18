#!/usr/local/bin/python
from fileCreation import *
from math import fabs
from primersPair import PrimersPair


class PseudoMultiplex:
    """
    This is a simple pseudo multiplex_machines who use targets saved in file (serialised) to not re-execute target.py script.

    :param number_of_targets: the number of target saved in the file (important for dill.load)
    """
    def __init__(self, number_of_targets):
        self.targets = self.sort_targets_by_bioimp(load_targets_from_file(number_of_targets))

    def create_chromosome(self):
        """
        Creates a chromosome with amplicon of the minimum size (TFGR minimum)

        :return: unique chromosome [primersPair0, primersPair1,...]\
        with primersPair having the lowest TFGR range for each target
        """

        chromosome = []
        for target in self.targets:
            chromosome.append(target.allPrimersPairs[0])
        return chromosome
    
    def run(self):
        """
        run the multiplex_machines

        :return: a chromosome with an optimize solution.
        """
        chromosome = self.check_condition()
        print(len(chromosome))
        return chromosome
    
    def check_condition(self):
        res = []
        for target in self.targets:
            add = True
            j = 0
            candidate = target.all_primers_pairs[j]
            if len(res) > 1:
                success = 0
                ok = True
                merged_at = -1
                while success < len(res) and ok:
                    for element in res:
                        if self.is_overlapping(candidate, element):
                            merge_primer_pair = self.merge_primerspairs(candidate, element)
                            if merge_primer_pair.TFGR[1] - merge_primer_pair.TFGR[0] <= 90:
                                merged_at = res.index(element)
                                candidate = merge_primer_pair
                                success += 1
                            else:
                                j += 1
                                if j < len(target.all_primers_pairs) - 1:
                                    candidate = target.all_primers_pairs[j]
                                    merged_at = -1
                                    success = 0
                                    break
                                else:
                                    ok = False
                                    add = False
                        else:
                            if self.check_primers_pairs(candidate, element):
                                # dimer detected or untargetted
                                j += 1
                                if j < len(target.all_primers_pairs) - 1:
                                    candidate = target.all_primers_pairs[j]
                                    merged_at = -1
                                    success = 0
                                    break
                                else:
                                    ok = False
                                    add = False
                            else:
                                success += 1
                                
                if merged_at != -1 and ok:
                    res[merged_at] = candidate
            if add:
                if candidate not in res:
                    res.append(candidate)
        return res

    def is_overlapping(self, primers_pair1, primers_pair2):
        """
        checks if two primers pairs are overlapping.

        :param primers_pair1: a primer pair object
        :param primers_pair2: a primer pair object

        :return: * True if the two primers pairs are overlapping
                 * False otherwise
        """
        if primers_pair1.left_primer.target.no_chromosome == primers_pair2.left_primer.target.no_chromosome:
            if (primers_pair2.TFGR[1] > primers_pair1.TFGR[1] > primers_pair2.TFGR[0]) and (
                            primers_pair1.TFGR[1] > primers_pair2.TFGR[0] > primers_pair1.TFGR[0]):
                return True
            elif (primers_pair2.TFGR[1] > primers_pair1.TFGR[0] > primers_pair2.TFGR[0]) and (
                            primers_pair1.TFGR[1] > primers_pair2.TFGR[1] > primers_pair1.TFGR[0]):
                return True
            elif primers_pair1.TFGR == primers_pair2.TFGR:
                return True
            else:
                return False
        else:
            return False

    def merge_primerspairs(self, primers_pair1, primers_pair2):
        """
        Merges two primers pairs when they overlap.

        :param primers_pair1: a primers pair object
        :param primers_pair2: a primers pair object

        :return: a new primer pairs object resulting of the addition of the two primers pairs.\
        i.e the left primer of the first and the right primer of the second.
        """
        primer_left = None
        primer_right = None
        if primers_pair1.TFGR[0] < primers_pair2.TFGR[0]:
            primer_left = primers_pair1.left_primer
        else:
            primer_left = primers_pair2.left_primer
        if primers_pair1.TFGR[1] > primers_pair2.TFGR[1]:
            primer_right = primers_pair1.right_primer
        else:
            primer_right = primers_pair2.right_primer
        merged_primers_pairs = PrimersPair(primer_left, primer_right)
        return merged_primers_pairs
    
    def compute_matches_between_primers(self, primer1, primer2):
        """
        Computes alignment between two primers.

        :param primer1: a Primer instance
        :param primer2: a Primer instance

        :return: the number of matches between the two primers
        """
        if len(primer1) > len(primer2):
            temp = primer2
            primer2 = primer1
            primer1 = temp
        bestMatches = 0
        mismatchedAssociated = 0
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        for i in range(1,len(primer1)+1):
            currentMatches = 0
            currentMismatched = 0
            p1temp = primer1[:i]
            p2temp = primer2[-i:]
            for j in range(i):
                if complement[p1temp[j]] == p2temp[j]:
                    # It's a match
                    currentMatches += 1
                else:
                    # It's a mismatch
                    currentMismatched += 1
            if currentMatches > bestMatches:
                bestMatches = currentMatches
                mismatchedAssociated = currentMismatched
        for i in range(len(primer1), 0, -1):
            currentMatches = 0
            currentMismatched = 0
            p1temp = primer1[-i:]
            p2temp = primer2[:i]
            for j in range(len(p1temp)):
                if complement[p1temp[j]] == p2temp[j]:
                    #It's a match
                    currentMatches += 1
                else:
                    #It's a mismatch
                    currentMismatched += 1
            if currentMatches > bestMatches:
                bestMatches = currentMatches
                mismatchedAssociated = currentMismatched
        #return bestMatches - mismatchedAssociated
        #Better ?
        return bestMatches
    
    def sort_targets_by_bioimp(self, targets):
        """
        Sorts targets by biological implication

        :param targets: a list containing Target instances

        :return: list of target sorted by important biological implication first.
        """
        damaging = []
        pot_damaging = []
        prob_poly = []
        unknown = []
        for target in targets:
            if target.bio_imp == "Damaging":
                damaging.append(target)
            elif target.bio_imp == "Probably Damaging":
                pot_damaging.append(target)
            elif target.bio_imp == "Probably Polymorphism":
                prob_poly.append(target)
            else:
                unknown.append(target)
        return damaging + pot_damaging + prob_poly + unknown

    def compute_score(self, primers_pair1, primers_pair2):
        """
        Computes the score of alignment of two primers (detection of dimers) and computes hybridization sites (true if it is, false otherwise)

        :param primers_pair1: a primersPair instance
        :param primers_pair2: a primersPair instance

        :return: * score_dim: a list containing the match scores (resulting of compute_matches_between_primers())\
                 for each possibility (4) for the two primers pairs. Score dimers

                 * check_amplification: a list containing boolean value resulting of search_hybridisation_sites function\
                 for each possibility between two primers pairs
        """
        scores_dim = [self.compute_matches_between_primers(primers_pair1.left_primer.sequence, primers_pair2.left_primer.sequence),
                      self.compute_matches_between_primers(primers_pair1.right_primer.sequence, primers_pair2.left_primer.sequence),
                      self.compute_matches_between_primers(primers_pair1.left_primer.sequence, primers_pair2.right_primer.sequence),
                      self.compute_matches_between_primers(primers_pair1.right_primer.sequence, primers_pair2.right_primer.sequence)]

        checks_amplification = [self.search_hybridisation_sites(primers_pair1.left_primer, primers_pair2.left_primer),
                                self.search_hybridisation_sites(primers_pair1.right_primer, primers_pair2.left_primer),
                                self.search_hybridisation_sites(primers_pair1.left_primer, primers_pair2.right_primer),
                                self.search_hybridisation_sites(primers_pair1.right_primer, primers_pair2.right_primer)]

        return scores_dim, checks_amplification

    def check_primers_pairs(self, primers_pair1, primers_pair2):
        """
        Checks if the new primers pairs is OK according to conditions

        :param primers_pair1: a primersPair instance
        :param primers_pair2: a primersPair instance

        :return: * True if the new primers pairs are compatible with the constraints

                 * False otherwise
        """
        scores = self.compute_score(primers_pair1, primers_pair2)
        if max(scores[0]) < 10 and True not in scores[1]:
            return False
        else:
            return True

    def search_hybridisation_sites(self, primer1, primer2):
        """
        Searches hybridization sites between two primers
        used by compute_score function

        :param primer1: a Primer instance
        :param primer2: a Primer instance

        :return: * True if there is hybridisation between two primers

                 * False otherwise
        """

        def check_hybridization(blast_result1, blast_result2):
            """
            Checks the hybridization of the two primers blast_result1 and blast_result2

            :param blast_result1: blast 'content' list
            :param blast_result2: blast 'content' list

            :return: * True if hybridisation

                     * False otherwise
            """
            if blast_result1[2] == 1 and blast_result2[2] == -1:
                if fabs(int(blast_result1[1]) - int(blast_result2[0])) < 1000:
                    return True
            elif blast_result1[2] == -1 and blast_result2[2] == 1:
                if fabs(int(blast_result2[1]) - int(blast_result1[0])) < 1000:
                    return True
            return False

        # check L-R and R-L
        for k1 in primer1.primer_hybridization_sites.keys():
            for k2 in primer2.primer_hybridization_sites.keys():
                if k1 == k2:
                    for blastResult1 in primer1.primer_hybridization_sites.get(k1):
                        for blastResult2 in primer2.primer_hybridization_sites.get(k2):
                            if check_hybridization(blastResult1, blastResult2):
                                return True
                    break
        return False
