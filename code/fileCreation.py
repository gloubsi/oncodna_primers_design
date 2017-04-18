#!/usr/local/bin/python

"""
File allowing to create all kind of useful files like saving targets's primers pair in bed file or
save target object into file for example.
"""
import dill
from config import *


def create_fasta_file(targets):
    """
    Creates a fasta file containing all sequences of targets

    :param targets: a list containing target objects.

    :return: nothing but creates a file sequences.fasta containing all targets's sequences.
    """
    with open(targets_sequences_file_path, "w") as f:
        for target in targets:
            f.write(">" + target.noChromosome + ":" + target.mutPosition + "\n" + target.sequence + "\n")


def create_fasta_file_primer(primer_seq):
    """
    Creates fasta file containing a sequence of primer

    :param primer_seq: a primer sequence

    :return: nothing but create primers.fasta containing the sequence of the primer passed in parameter.
    """
    with open(primer_to_blast_file_path, "w") as f:
            f.write(">primer\n" + primer_seq + "\n")


def create_bed_file(targets):
    """
    Creates a bed file containing chromosomal range based on the first primers pair of each target.
    chromosomal range: last position of the left primer and first position of the right

    :param targets: a list containing target objects.

    :return: nothing but create chrom.bed a file containing such chromosomal range
    """
    c = 1
    with open(chrom_file_path, "w") as b:
        for target in targets:
            if target.all_primers_pairs:
                best_primers_pair = target.all_primers_pairs[0]
                b.write("#" + str(c) + "\n#bioImp: " + str(target.bio_imp) +
                        "\n#varID: " + target.var_key +
                        "\n" + best_primers_pair.no_chromosome + "\t")
                for position in best_primers_pair.TFGR:
                    b.write(str(position) + "\t")
                b.write("\n")
                c += 1


def create_bed_file_chromosome(chromosome):
    """
    Creates a bed file with a list a primers pairs.
    Originally created for the multiplex_machines algorithm.

    :param chromosome: a list containing primers pair object.

    :return: nothing but create a file chrom.bed containing chromosomic range according to primers pairs.
    """
    c = 1
    with open(chrom_file_path, "w") as b:
        for primersPair in chromosome:
            b.write("#" + str(c) + "\n" + "#mutation position: " + str(primersPair.left_primer.target.mutation_pos) + "\n")
            b.write(primersPair.left_primer.target.no_chromosome + "\t" + str(primersPair.TFGR[0]) + "\t" + str(primersPair.TFGR[1]) + "\n")
            c += 1


def save_targets_into_file(targets):
    """
    Saves targets into files. Important for not generating again primers pairs of targets.

    :param targets: a list containing target object

    :return: nothing but create targets_datas.pkl containing target objects (serialisation)
    """
    with open(serialisation_file_path, "wb") as output:
        for target in targets:
            dill.dump(target, output)


def load_targets_from_file(numberOfTargets):
    """
    Loads targets from the file targets_datas.pkl

    :param numberOfTargets: the number of targets stores in the file

    :return: a list of the target objects that have been saved previously
    """
    targets = []
    with open(serialisation_file_path, 'rb') as f:
        for _ in range(numberOfTargets):
            targets.append(dill.load(f))
    return targets
