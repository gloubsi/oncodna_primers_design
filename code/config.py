# File that need to be created before running the script
genome_file_path = "/home/bioinfo/opt/reference/hg19/hg19_chrOnly.fa"
# file path creation:
serialisation_file_path = "/home/npotie/stage/datas/targets_datas.pkl"
chrom_file_path = "/home/npotie/stage/datas/chrom.bed"
primer_to_blast_file_path = "/home/npotie/stage/datas/primers.fasta"
blast_output_file_path = "/home/npotie/stage/datas/outputBlast.txt"
targets_sequences_file_path = "/home/npotie/stage/datas/sequences.fasta"

# Launch (from main.py) parameters. Do not make multiplex and pseudo_multiplex both True
# if overlap = True, no multiplex algorithms will be run.
blast = True
overlap = True

# multiplex parameters:
multiplex = False
size_of_pop = 100
nb_gen = 500
# mutation probability
PM = 0.2
# Cross-over probability
PC = 0.3

# pseudo_multiplex parameters:
pseudo_multiplex = False

# Target parameters.
# Don't make sort_tfgr and sort_by_penalty both True.
seq_range = 100
sort_tfgr = True
sort_by_penalty = False
remove_hybridisation = False

# Primer parameters
gmaf_limit = 0.05

# Primer3 parameters:
PRIMER_OUTSIDE_PENALTY = 5

PRIMER_NUM_RETURN = 10
PRIMER_PAIR_NUM_RETURNED = 10

PRIMER_OPT_SIZE = 20
PRIMER_MIN_SIZE = 18
PRIMER_MAX_SIZE = 25

PRIMER_OPT_TM = 60.0
PRIMER_MIN_TM = 57.0
PRIMER_MAX_TM = 63.0

PRIMER_MIN_GC = 30.0
PRIMER_MAX_GC = 70.0

PRIMER_SALT_MONOVALENT = 50.0
PRIMER_DNA_CONC = 50.

PRIMER_MAX_POLY_X = 5
PRIMER_MAX_SELF_ANY = 4
PRIMER_MAX_NS_ACCEPTED = 0
PRIMER_MAX_SELF_END = 3

PRIMER_PAIR_MAX_COMPL_END = 3
PRIMER_PAIR_MAX_COMPL_ANY = 4

PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1
PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = 1

PRIMER_PRODUCT_SIZE_RANGE = [[51, 140]]

# BLAST parameters
blast_executable = "/home/npotie/blast/bin/blastn"
e_value_limit = 1
perc_identity = 100
p_align = 80

# Data file for tests
data_1 = "/home/npotie/stage/datas/data_1.csv"
data_16 = "/home/npotie/stage/datas/data_16.csv"
data_21 = "/home/npotie/stage/datas/data_21.csv"
data1_10 = "/home/npotie/stage/datas/data1_10.csv"