#!/usr/local/bin/python

import random


class Chromosome:
    """
    A chromosome is one 'entity' of the population of the genetic algorithm

    :param components_list: a list containing primersPair instance (size equal to the number of targets)
    """

    def __init__(self, components_list):
        self.components = components_list
        self.fitness = self.compute_chromosome_fitness()
    
    def __str__(self):
        return str(self.components) + "\nfitness: " + str(self.fitness)
    
    def __eq__(self, other):
        return self.components == other.components
    
    def mutate(self, no_target, new_primers_pair):
        """
        Mutate the chromosome, change the primer pair of the noTarget by the new primers pair.

        :param no_target: the target number
        :param new_primers_pair: the new instance of primersPair at the index no_target

        :return: nothing but have a side effect on components at index no_target.
        """
        self.components[no_target] = new_primers_pair
    
    def compute_chromosome_fitness(self):
        """
        Computes the fitness of the whole chromosome.

        :return: a fitness, wich is a float value\
        corresponding to the addition of the fitness all two primers pairs of the components list.
        """
        fitness = 0
        for i in range(len(self.components)-1):
            primers_pair1 = self.components[i]
            for j in range(i+1, len(self.components)):
                primers_pair2 = self.components[j]
                fitness += self.compute_fitness_between_primers(primers_pair1.left_primer.sequence, primers_pair2.left_primer.sequence)
                fitness += self.compute_fitness_between_primers(primers_pair1.left_primer.sequence, primers_pair2.right_primer.sequence)
                fitness += self.compute_fitness_between_primers(primers_pair1.right_primer.sequence, primers_pair2.left_primer.sequence)
                fitness += self.compute_fitness_between_primers(primers_pair1.right_primer.sequence, primers_pair2.right_primer.sequence)
        return fitness
    
    def update_fitness(self):
        """
        Updates the fitness of a chromosome

        :return:  nothing by have a side effect on fitness of the chromosome
        """
        self.fitness = self.compute_chromosome_fitness()
    
    def cross_over(self, other):
        """
        Performs a crossing-over between two chromosomes: self and other.
        if the cutting point < 12 (half of length), exchange components from the front end part of the two chromosomes
        if cutting point > 12, exchange components from the rear end part of the two chromosomes
        So if cutting point = 2, exchange the first 2 components. If cutting point = 22, exchange the three last components.

        :param other: a other instance of chromosome object

        :return: nothing but performs the cross-over
        """
        cutting_point = random.randint(1, 24)
        if cutting_point <= 12:
            temp = self.components[:cutting_point]
            self.components[:cutting_point] = other.components[:cutting_point]
            other.components[:cutting_point] = temp
        else:
            temp = self.components[-(len(self.components) - cutting_point):]
            self.components[-(len(self.components) - cutting_point):] = other.components[-(len(self.components) - cutting_point):]
            other.components[-(len(self.components) - cutting_point):] = temp

    def compute_fitness_between_primers(self, primer1, primer2):
        """
        Computes alignment between two primers.

        :param primer1: a Primer instance
        :param primer2: a Primer instance

        :return: the number of mismatched - matches
        """
        if len(primer1) > len(primer2):
            temp = primer2
            primer2 = primer1
            primer1 = temp
        best_matches = 0
        mismatched_associated = 0
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        for i in range(1,len(primer1)+1):
            current_matches = 0
            current_mismatched = 0
            p1temp = primer1[:i]
            p2temp = primer2[-i:]
            for j in range(i):
                if complement[p1temp[j]] == p2temp[j]:
                    # It's a match
                    current_matches += 1
                else:
                    # It's a mismatch
                    current_mismatched += 1
            if current_matches > best_matches:
                best_matches = current_matches
                mismatched_associated = current_mismatched
        for i in range(len(primer1), 0, -1):
            current_matches = 0
            current_mismatched = 0
            p1temp = primer1[-i:]
            p2temp = primer2[:i]
            for j in range(len(p1temp)):
                if complement[p1temp[j]] == p2temp[j]:
                    # It's a match
                    current_matches += 1
                else:
                    # It's a mismatch
                    current_mismatched += 1
            if current_matches > best_matches:
                best_matches = current_matches
                mismatched_associated = current_mismatched
        # return best_matches - mismatched_associated
        # Better ?
        return mismatched_associated - best_matches
