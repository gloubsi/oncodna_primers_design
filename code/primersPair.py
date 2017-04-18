#!/usr/local/bin/python


class PrimersPair:
    """
    PrimersPair is an object that allow to link a left primer and a right primer that amplify a region
    (containing a region)

    Two primers pairs are said equal if their left primers and right primers are equal respectively

    :param left_primer: the left primer instance
    :param right_primer: the right primer instance
    """
    def __init__(self, left_primer, right_primer):
        self.left_primer = left_primer
        self.right_primer = right_primer
        self.TFGR = [self.left_primer.TFGP, self.right_primer.TFGP]
        self.penalty = left_primer.penalty + right_primer.penalty
        # amplicon_range between primers
        self.amplicon_range = self.TFGR[1] - self.TFGR[0]
    
    def __str__(self):
        return "\n\nposition of mutation: " + str(self.left_primer.target.mutation_pos) + \
               "\nleft primer: " + str(self.left_primer) + \
               "\nright primer: " + str(self.right_primer) + \
               "\namplicon size: " + str(self.amplicon_range)
    
    def __eq__(self, other):
        """
        Two primers pairs are said equal if their left primers and right primers are equal respectively

        :param other: an another primers pair instance

        :return: * True if the two primers pair are equal

                 * False otherwise
        """
        return (self.left_primer == other.left_primer) and (self.right_primer == other.right_primer)
        
    def __ne__(self, other):
        return not self.__eq__(other)
