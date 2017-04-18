#!/usr/local/bin/python
from fileCreation import *
from parsing import get_data_from_csv
from postTreatments import PostTreatments
from config import *
import sys
import time
from functools import reduce


def launch(path_of_data_file):
    """
    Launches the algorithm.

    :param path_of_data_file: the path of the data file

    :return: nothing but creates a bed file corresponding to good chromosomic range.
    """
    print("Welcome to primers design program!")
    print("-------------------------------------------------------------")
    print("Documentation of the program in pydoc/_build/html/index.html")
    print("-------------------------------------------------------------")
    print("Processing...")
    start_time = time.time()
    targets = get_data_from_csv(path_of_data_file)
    size = len(targets)
    res = []
    good_targets = []
    for target in targets:
        target.compute_list_primers_pairs()
        if target.all_primers_pairs:
            print("target design successfully: " + target.no_chromosome + " " + str(target.mutation_pos) + " " + target.bio_imp)
            res.append(target.all_primers_pairs[0])
            good_targets.append(target)
        else:
            print("fail to design target: " + target.no_chromosome + " " + str(target.mutation_pos) + " " + target.bio_imp)
    if len(good_targets) > 1:
        pt = PostTreatments(good_targets)
        if overlap:
            print("checking overlap ...")
            res = pt.check_overlap()
        elif size > 1 and multiplex:
            print("multiplex in progress ...")
            res = pt.run_multiplex()
        elif size > 1 and pseudo_multiplex:
            print("pseudo multiplex in progress ...")
            res = pt.run_pseudo_multiplex()

    if not multiplex:
        create_bed_file_chromosome(res)
    if multiplex:
        create_bed_file_chromosome(res.components)
    print("program finish in: " + seconds_to_str(time.time() - start_time))


def seconds_to_str(t):
    """
    Converts the time from seconds to hh:mm:ss:ms
    :param t: the time in second
    :return: return a string as hh:mm:ss:ms
    """
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(t*1000,), 1000, 60, 60])


if __name__ == '__main__':
    launch(sys.argv[1])
