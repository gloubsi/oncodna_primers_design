from primersPair import *
from multiplex_machines.multiplex import MultiplexMachine
from fileCreation import save_targets_into_file
from multiplex_machines.pseudoMultiplex import PseudoMultiplex
from config import *


class PostTreatments:
    """
    Class that make post treatments i.e checking and merge overlapping best primers pairs of 2 distinct targets,\
    the multiplex_machines algorithm...
    All operations that we can do when all targets have been generated.

    **We suppose that the best primers pair for each target is the first in the list**

    :param targets: a list containing target instance (when all_primers_pairs has been computed)

    :ivar l_targets: a list of Target
    """
    def __init__(self, targets):
        save_targets_into_file(targets)
        self.l_targets = targets

    def check_overlap(self):
        """
        When two targets are overlapping each other i.e if his first primers pair are overlapping
        then merge these primers pairs.

        :return: a list containing the first primers pair of each target, no overlap.\
        so one primers pair car refer to one or more targets.
        """
        res = [self.l_targets[0].all_primers_pairs[0]]
        for target in self.l_targets[1:]:
            overlap = False
            same_chr = True
            for element in res:
                if target.no_chromosome == element.right_primer.target.no_chromosome:
                    same_chr = True
                    if self.is_overlapping(target.all_primers_pairs[0], element):
                        overlap = True
                        merge_target = self.merge_primerspairs(target.all_primers_pairs[0], element)
                        if merge_target.TFGR[1] - merge_target.TFGR[0] <= 90:
                            # print("merged ! " + merge_target.left_primer.target.no_chromosome)
                            res[res.index(element)] = merge_target
                            break
                else:
                    same_chr = False

            if not overlap:
                # print("no overlapping: " + target.no_chromosome)
                res.append(target.all_primers_pairs[0])
            elif not same_chr:
                # print("not same chr: " + target.no_chromosome)
                res.append(target.all_primers_pairs[0])
        return res

    def is_overlapping(self, primers_pair1, primers_pair2):
        """
        checks if two primers pairs are overlapping.

        :param primers_pair1: a primer pair object
        :param primers_pair2: a primer pair object

        :return: * True if the two primers pairs are overlapping

                 * False otherwise
        """
        if (primers_pair2.TFGR[1] >= primers_pair1.TFGR[1] >= primers_pair2.TFGR[0]) and (
                primers_pair1.TFGR[1] >= primers_pair2.TFGR[0] >= primers_pair1.TFGR[0]):
            return True
        elif (primers_pair2.TFGR[1] >= primers_pair1.TFGR[0] >= primers_pair2.TFGR[0]) and (
                primers_pair1.TFGR[1] >= primers_pair2.TFGR[1] >= primers_pair1.TFGR[0]):
            return True
        elif primers_pair1.TFGR == primers_pair2.TFGR:
            return True
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
        if primers_pair1.TFGR[0] <= primers_pair2.TFGR[0]:
            primer_left = primers_pair1.left_primer
        else:
            primer_left = primers_pair2.left_primer
        if primers_pair1.TFGR[1] >= primers_pair2.TFGR[1]:
            primer_right = primers_pair1.right_primer
        else:
            primer_right = primers_pair2.right_primer
        merged_primers_pairs = PrimersPair(primer_left, primer_right)
        return merged_primers_pairs

    def run_multiplex(self):
        """
        Runs the multiplex_machines algorithm
        Be careful this algorithm are not really optimized and delete too many targets !

        :return: a chromosome of optimized primers_pair
        """
        multiplex = MultiplexMachine(len(self.l_targets))
        return multiplex.run()

    def run_pseudo_multiplex(self):
        """
        Runs the pseudo multiplex_machines algorithm
        Be careful this algorithm are not really optimized and delete too many targets !

        :return: a chromosome of optimized primers_pair
        """
        p_multiplex = PseudoMultiplex(len(self.l_targets))
        return p_multiplex.run()
