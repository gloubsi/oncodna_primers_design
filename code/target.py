#!/usr/local/bin/python

import primer3
import pysam
import myvariant
import math

from primer import Primer
from primersPair import PrimersPair
from blast import *
from config import *


class Target:
    """
    Target is an object that refer to one line in the CSV that is to say a variant.

    :param target_components: a list resulting of the command line.split() in the csv. Contains:

    :ivar no_chromosome: the chromosome number where the variant is.

    :ivar mutation_pos: the position of the mutation (variant).

    :ivar var_freq: the frequency of the variant.

    :ivar bio_imp: the biological implication of the variant could be (order by priority).

                   * Damaging

                   * Potentially damaging

                   * Probably polymorphism

                   * Unknown

    :ivar var_key: the id of the variant containing chromosome number + mutation position.

    :ivar range: the number of nucleotides taken into account around variant\
                 for designing primers pairs.

    :ivar sequence: the target sequence containing the variant (empty initially)

    :ivar all_primers_pairs: a list containing all valid primers pairs within the target sequence\
                             empty at initialisation
    """
    def __init__(self, target_components):
        self.no_chromosome = target_components[3]
        self.mutation_pos = int(target_components[4])
        self.var_freq = float(target_components[9])
        self.bio_imp = target_components[11]
        self.var_key = target_components[16]
        self.range = seq_range
        self.sequence = ""
        self.all_primers_pairs = []
            
    def __str__(self):
        """
        Fancy printing the target object with all her useful proprieties
        """
        return "Sequence: " + self.sequence + \
               "\nChromosome number: " + self.no_chromosome + \
               "\nMutation Position: " + str(self.mutation_pos) + \
               "\nBiological impact: " + self.bio_imp + \
               "\nNumber of primers pairs: " + str(len(self.all_primers_pairs))
    
    def get_target_sequence(self, range):
        """
        Fetches from the fasta file hg19 (saved locally) the sequence containing the variant.
        pysam library is used.

        :param range: the number of nucleotide taken into account around the variant

        :return: the sequence containing the variant (in the middle of the sequence)
        """
        fasta_file = pysam.Fastafile(genome_file_path)
        sequence = fasta_file.fetch(self.no_chromosome, self.mutation_pos - range, self.mutation_pos + range + 1)
        self.sequence = sequence

    def design_primers(self, range, sequence):
        """
        Designs a pair of primers (right and left) using primer3 library.
        Note that primer3 won't design primers with hairpin, self hybridisation and
        hybridisation with the right and left primers
        Be careful, parameters are approximated to be as correct as Thermofisher.
        See config.py for parametrisation.

        :param range: the number of nucleotide taken into account around the variant

        :param sequence: the sequence containing the variant obtained by get_target_sequence function.

        :return: a primers pair
        """
        primer = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_TEMPLATE': sequence
                    },
                    {
                        # type of procedure. Get the primers between the SNP
                        'PRIMER_TASK': 'pick_sequencing_primers',
                        'SEQUENCE_TARGET': [range + 1, 1],
                        'PRIMER_OUTSIDE_PENALTY': PRIMER_OUTSIDE_PENALTY,
                        
                        # Maximum number of pairs of primers returned
                        # Default value: 5
                        'PRIMER_NUM_RETURN': PRIMER_NUM_RETURN,
                        'PRIMER_PAIR_NUM_RETURNED': PRIMER_PAIR_NUM_RETURNED,
                        
                        # size of the primer
                        'PRIMER_OPT_SIZE': PRIMER_OPT_SIZE,
                        'PRIMER_MIN_SIZE': PRIMER_MIN_SIZE,
                        'PRIMER_MAX_SIZE': PRIMER_MAX_SIZE,
                        
                        # Temperature of the primer
                        'PRIMER_OPT_TM': PRIMER_OPT_TM,
                        'PRIMER_MIN_TM': PRIMER_MIN_TM,
                        'PRIMER_MAX_TM': PRIMER_MAX_TM,
                        
                        # GC propensity
                        'PRIMER_MIN_GC': PRIMER_MIN_GC,
                        'PRIMER_MAX_GC': PRIMER_MAX_GC,
                        
                        # The maximum allowable length of a mononucleotide repeat (AAAAA for example)
                        # Default value: 5
                        'PRIMER_MAX_POLY_X': PRIMER_MAX_POLY_X,
                        
                        # The millimolar (mM) concentration of monovalent salt cations. Used to calculate the melting T.
                        # Default value: 50
                        'PRIMER_SALT_MONOVALENT': PRIMER_SALT_MONOVALENT,
                        
                        # A value to use as nanomolar (nM) concentration of each annealing oligo over the course the PCR
                        # Primer3 uses this argument to esimate oligo melting temperatures.
                        # Default value: 50
                        'PRIMER_DNA_CONC': PRIMER_DNA_CONC,
                        
                        # Maximum number of unknown bases (N) allowable in any primer
                        # Default value: 0
                        'PRIMER_MAX_NS_ACCEPTED': PRIMER_MAX_NS_ACCEPTED,
                        
                        # Maximum allowable 3'-anchored global alignment score
                        # when testing a single primer for self-complementarity.
                        # Default value: 3
                        'PRIMER_MAX_SELF_END': PRIMER_MAX_SELF_END,
                        
                        # tries to bind the 3-END of the left primer to the right primer
                        # and scores the best binding it can find.
                        # Similar to PRIMER_MAX_SELF_END
                        # Default value: 3
                        'PRIMER_PAIR_MAX_COMPL_END': PRIMER_PAIR_MAX_COMPL_END,
                        
                        # primer3 will use thermodynamic models to calculate
                        # the propensity of oligos to form hairpins and dimers.
                        # Check for hairpin et dimerization between primers pairs.
                        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT,
                        
                        # primer3 will use thermodynamic models to calculate the the propensity of oligos
                        # to anneal to undesired sites in the template sequence.
                        'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT,
                        
                        # describes the tendency of a primer to bind to itself. Check for self-dimerization
                        'PRIMER_MAX_SELF_ANY': PRIMER_MAX_SELF_ANY,
                        
                        # Describes the tendency of the left primer to bind to the right primer.
                        # It is the maximum allowable local alignment score when testing for complementarity
                        # between left and right primers. Similar to PRIMER_MAX_SELF_ANY.
                        # Default value: 8
                        'PRIMER_PAIR_MAX_COMPL_ANY': PRIMER_PAIR_MAX_COMPL_ANY,
                        
                        # Require the specified number of consecutive
                        # Gs and Cs at the 3' end of both the left and right primer.
                        # (This parameter has no effect on the internal oligo if one is requested.)
                        # 'PRIMER_GC_CLAMP': 0,
                        
                        # The maximum allowable length of a mononucleotide repeat, for example AAAAAA.
                        # Don't impact the success
                        # 'PRIMER_MAX_POLY_X': 5,
                        
                        # Length PCR product result. We want a product between 50 and 100 base pairs.
                        # Default value: 100-300
                        'PRIMER_PRODUCT_SIZE_RANGE': PRIMER_PRODUCT_SIZE_RANGE,
                    })
                )
        return primer
    
    def compute_list_primers_pairs(self):
        """
        Main function to design a set of primers pairs that amplify the variant.
        Remind you that there is range nucleotide around the variant.
        The function does the following:

        #. compute_optimal_primers_pairs() computes all possible primers pairs of the sequence.\
        to do that, each iteration reduce the length of the  (non symmetrical sub-sequence calculation) by 5\
        and design_primers() is execute at each iteration.\
        Note that we keep just primers pairs that produce an amplicon of length <= 90 (requirement from Thermofisher),\
        and no snp in the primers sequences.\
        We will get a list containing all primers pairs calculated.

        #. perform_blast() performs a blast for each primers pair and each primer taken individually.\
        will add hybridisation site for each primer.\
        Refers to blast class for more information.

        #. sort_by_tfgr() a function that sort primers pairs by size of amplicon (minimum range of tfgr).

        #. remove_primers_pair_with_hybridisation() will remove all primers pairs whose primers\
        could hybridise something elsewhere in the genome.\
        Be careful, it's BLAST class that initialize hybridisation sites

        :param blast: * True if the user want to blast with all primers pairs computed (default)

                      * False otherwise

        refer to config.py to set boolean

        :return: nothing but store the result in all_primers_pairs variable of target.
        """

        self.all_primers_pairs = self.compute_optimal_primers_pairs()
        if blast:
            self.perform_blast()
        if sort_tfgr:
            self.sort_by_tfgr()
        if sort_by_penalty:
            self.sort_by_penalty()
        if remove_hybridisation:
            self.remove_primers_pair_with_hybridisation()

    def compute_optimal_raw_primers_pairs(self):
        """
        /!\ do not call this function to compute all_primers_pairs final list
        call self.compute_optimal_primers_pairs() for no post-treatment
        or call self.compute_list_primers_pairs() for some post-treatment.

        Computes list of primers pairs within the sequence containing the variant.
        Be careful primers pair are dic (raw output from primer3 algo) and not
        primersPair instance.

        :return all_primers_pairs: a list containing raw output from design_primers function such as :\
                                    [primers1, primers2,...]\
                                    with primers 1 a dic (output of primer3 algo)\
        :return ranges: containing range_temp_left for each primers added in all_primers_pairs.
        """
        all_primers_pairs = []
        range_temp_left = self.range
        finish = False
        self.get_target_sequence(self.range)
        ranges_left = []
        while not finish:
            while range_temp_left > 25:
                range_temp_left -= 5
                range_temp_right = self.range
                while range_temp_right > 25:
                    range_temp_right -= 5
                    sequence = self.sequence[(self.range - range_temp_left):-(self.range - range_temp_right)]
                    primers = self.design_primers(range_temp_left, sequence)
                    if "PRIMER_LEFT_0_SEQUENCE" in primers and "PRIMER_RIGHT_0_SEQUENCE" in primers:
                        variant_position = range_temp_left + 1
                        if primers["PRIMER_LEFT_0"][0] + primers["PRIMER_LEFT_0"][1] < variant_position < \
                                primers["PRIMER_RIGHT_0"][0] - primers["PRIMER_RIGHT_0"][1]:

                            # sequence_bt_primers = sequence[primers["PRIMER_LEFT_0"][0] + primers["PRIMER_LEFT_0"][1]: primers["PRIMER_RIGHT_0"][0] - primers["PRIMER_RIGHT_0"][1]]
                            all_primers_pairs.append(primers)
                            ranges_left.append(range_temp_left)
            finish = True
        return all_primers_pairs, ranges_left

    def compute_optimal_primers_pairs(self):
        """
        Computes a list of primers pairs containing primersPair instance
        Use this function to get a list of primers pairs respecting all constraints
        with no post-treatment (no BLAST, no sort,...)

        :return: a list of primersPair instance.
        """
        all_primers_pairs = []
        result = self.compute_optimal_raw_primers_pairs()
        all_raw_primers_pairs = result[0]
        ranges_left = result[1]
        SNP = self.get_snp()
        for i in range(len(all_raw_primers_pairs)):

            left_primer = Primer(all_raw_primers_pairs[i]["PRIMER_LEFT_0_SEQUENCE"], "left",
                                 all_raw_primers_pairs[i]["PRIMER_LEFT_0"][0],
                                 ranges_left[i],
                                 self,
                                 all_raw_primers_pairs[i]["PRIMER_LEFT_0_PENALTY"])

            right_primer = Primer(all_raw_primers_pairs[i]["PRIMER_RIGHT_0_SEQUENCE"], "right",
                                  all_raw_primers_pairs[i]["PRIMER_RIGHT_0"][0],
                                  ranges_left[i],
                                  self,
                                  all_raw_primers_pairs[i]["PRIMER_RIGHT_0_PENALTY"])

            if not left_primer.dinucleotides_repeats() and not right_primer.dinucleotides_repeats():
                if math.fabs(right_primer.TFGP - left_primer.TFGP) <= 90:
                    if right_primer.check_snp_in_primer(SNP) == 0 and \
                            left_primer.check_snp_in_primer(SNP) == 0:

                        primers_pair = PrimersPair(left_primer, right_primer)
                        if primers_pair not in all_primers_pairs:
                            all_primers_pairs.append(primers_pair)

        self.all_primers_pairs = all_primers_pairs
        return all_primers_pairs

    def get_snp(self):
        """
        Gets all existing snp (from my variant library) in the target sequence.
        Why ? Because a snp can't be in the primers sequences.
        internet connection required.

        :return: a list of snp such as [[id, position, gmaf],[...],[...],...]
        """
        snp_info = []
        mv = myvariant.MyVariantInfo()
        res = mv.query(self.no_chromosome + ":" + str(int(self.mutation_pos) - self.range) +
                       "-" + str(int(self.mutation_pos) + self.range), fields='dbsnp', size=1000)
        for element in res["hits"]:
            if "dbsnp" in element:
                if "gmaf" in element["dbsnp"]:
                    snp_info.append([element["dbsnp"]["rsid"], element["dbsnp"]["hg19"]["start"], element["dbsnp"]["gmaf"]])
        return snp_info

    def perform_blast(self):
        """
        Performs a blast for each primers pair and each primer taken individually\
        will add hybridisation site for each primer.\
        Refers to blast class for more information.\
        The attribute all_primers_pairs must **not** be empty (run compute_list_primers_pairs() before)

        :return: nothing but have a side effect on all_primers_pairs
        """
        if self.all_primers_pairs:
            b = Blast(self.all_primers_pairs)
            b.add_hybridization_sites()
            self.all_primers_pairs = b.all_primers_pairs

    def sort_by_tfgr(self):
        """
        Sorts primers pairs by size of amplicon (minimum range of tfgr)\
        The attribute all_primers_pairs must **not** be empty (run compute_list_primers_pairs() before).
        Bubble sort.

        :return: nothing but have a side effect on all_primers_pairs.
        """
        if self.all_primers_pairs:
            for i in range(len(self.all_primers_pairs)):
                for j in range(len(self.all_primers_pairs) - 1 - i):
                    if self.all_primers_pairs[j].amplicon_range > self.all_primers_pairs[j + 1].amplicon_range:
                        self.all_primers_pairs[j], self.all_primers_pairs[j + 1] = self.all_primers_pairs[j + 1], self.all_primers_pairs[j]
                    elif self.all_primers_pairs[j].amplicon_range == self.all_primers_pairs[j + 1].amplicon_range:
                        if self.all_primers_pairs[j].TFGR[1] > self.all_primers_pairs[j + 1].TFGR[1]:
                            self.all_primers_pairs[j], self.all_primers_pairs[j + 1] = self.all_primers_pairs[j + 1], self.all_primers_pairs[j]

    def sort_by_penalty(self):
        """
        Sorts primers pairs by ascending penalties.
        Primer3 compute for each left primer and right primer a penalty.
        Higher is the penalty further the characteristics from optimal parameters passed in primer3 algorithm.
        Bubble sort

        :return: nothing but have a side effect on all_primers_pairs.
        """
        if self.all_primers_pairs:
            for i in range(len(self.all_primers_pairs)):
                for j in range(len(self.all_primers_pairs) - 1 - i):
                    if self.all_primers_pairs[j].penalty >= self.all_primers_pairs[j+1].penalty:
                        self.all_primers_pairs[j], self.all_primers_pairs[j + 1] = self.all_primers_pairs[j + 1], self.all_primers_pairs[j]

    def remove_primers_pair_with_hybridisation(self):
        """
        Removes all primers pairs whose primers could hybridise something elsewhere in the genome.\
        The attribute all_primers_pairs must **not** be empty (run compute_list_primers_pairs() before)
        Be careful, it's BLAST class that initialize hybridisation sites

        :return: nothing but have a side effect on all_primers_pairs
        """
        if self.all_primers_pairs:
            new_primers_pairs = []
            for primers_pair in self.all_primers_pairs:
                if not primers_pair.right_primer.primer_hybridization_sites and not primers_pair.left_primer.primer_hybridization_sites:
                    new_primers_pairs.append(primers_pair)
            self.all_primers_pairs = new_primers_pairs


