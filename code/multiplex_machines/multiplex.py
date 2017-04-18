#!/usr/local/bin/python
import random
from itertools import groupby
from operator import attrgetter

from fileCreation import *
from multiplex_machines.chromosome import Chromosome
from config import *


class MultiplexMachine:
    """
    multiplex_machines machine using genetic algorithm.
    !!! side effect on targets list

    :param number_of_targets: the number of target saved in the file (important for dill.load)
    """
    def __init__(self, number_of_targets):
        self.targets = load_targets_from_file(number_of_targets)

    def run(self):
        """
        Runs the multiplex_machines algorithm

        :return: a chromosome object whose components are optimised with our constraints (best fitness value)
        """
        population = self.create_initial_pop()
        print("Population generated!")
        current_generation = 0
        # number of chromosomes selected by roulette wheel pooled in the new population.
        nb_best_chromosomes = []
        # Checks occurences of the number of best chromosomes repetition
        occurences_b_c = []
        # Gets just the repetition
        occurences = []
        # break if we exceed number of generation or 5 times the same number of best chromosomes.
        while (current_generation < nb_gen) and (len(occurences) == 0):
            best_chromosomes = []
            occurences_b_c = [(k, sum(1 for _ in g)) for k, g in groupby(nb_best_chromosomes)]
            occurences = [x[1] for x in occurences_b_c if x[1] > 5]
            # print(occurences_b_c)
            # print(occurences)
            for i in range(len(population)):
                is_mutation = False
                is_cross_over = False
                chromosomes_selected = self.roulette_wheel_selection(population)
                x_prime = Chromosome(chromosomes_selected[0].components)
                y_prime = Chromosome(chromosomes_selected[1].components)
                rPm = random.random()
                if rPm <= PM:
                    self.perform_mutation(x_prime)
                    self.perform_mutation(y_prime)
                    is_mutation = True
                rPc = random.random()
                if rPc <= PC:
                    x_prime.cross_over(y_prime)
                    is_cross_over = True
                if is_mutation or is_cross_over:
                    x_prime.update_fitness()
                    y_prime.update_fitness()
                    if x_prime not in best_chromosomes:
                        best_chromosomes.append(x_prime)
                    if y_prime not in best_chromosomes:
                        best_chromosomes.append(y_prime)
                    
            nb_best_chromosomes.append(len(best_chromosomes))
            new_population =[x for x in best_chromosomes]
            while len(new_population) < size_of_pop:
                best_fitness_chr = max(population, key=attrgetter('fitness'))
                new_population.append(best_fitness_chr)
                population.remove(best_fitness_chr)
            population = new_population
            current_generation += 1
        
        return max(population, key=attrgetter('fitness'))
        
    def create_initial_pop(self):
        """
        Creates chromosomes into population for the multiplex_machines genetic algorithm.
        A chromosome is composed of targets (sequences) T that is (G, Pf, Pr, var_key).
        The melting temperature of primers for the same target are roughly the same (done by primer3).

        :return: chromosome containing (G, Pf, Pr): G is the pool number (always 1 is this case),\
        Pf the number of left primers, Pr the number of right primers
        """
        population =[]
        for _ in range(size_of_pop - 1):
            chromosome = self.create_chromosome()
            if chromosome not in population:
                population.append(chromosome)
        population.append(self.best_chromosome_simple_plex())
        return population
    
    def best_chromosome_simple_plex(self):
        """
        :return: chromosome containing primers pairs closest to the variant in each target
        """
        chromosome_components = []
        for i in range(len(self.targets)):
            chromosome_components.append(self.targets[i].all_primers_pairs[0])
        return Chromosome(chromosome_components)
    
    def create_chromosome(self):
        """
        Generates one chromosome object with random primers pairs as components

        :return: a chromosome containing primers pairs.
        """
        chromosome_components = []
        for i in range(len(self.targets)):
            randomPrimersPair = random.choice(self.targets[i].all_primers_pairs)
            chromosome_components.append(randomPrimersPair)
        return Chromosome(chromosome_components)
    
    def perform_mutation(self, chromosome):
        """
        Performs a mutation. Pick a random target, pick a random primer pair of the given target.
        The new primers pair picked will replaced the primers pair in place

        :param chromosome: a chromosome object instance

        :return: the chromosome with the random mutation applied
        """
        new_chromosome = Chromosome(chromosome.components)
        random_no_target = random.randint(0, len(self.targets) - 1)
        random_primers_pair = random.choice(self.targets[random_no_target].all_primers_pairs)
        new_chromosome.mutate(random_no_target, random_primers_pair)
        while new_chromosome.fitness < chromosome.fitness:
            random_no_target = random.randint(0, len(self.targets) - 1)
            random_primers_pair = random.choice(self.targets[random_no_target].all_primers_pairs)
            new_chromosome = Chromosome(chromosome.components)
            new_chromosome.mutate(random_no_target, random_primers_pair)
        chromosome.mutate(random_no_target, random_primers_pair)
    
    def roulette_wheel_selection(self, population):
        """
        Selects two chromosomes by the roulette wheel selection method.

        :param population: a list containing all chromosomes instances of the population.

        :return: two chromosomes selected by roulette wheel selection method
        """
        global_fitness = 0.
        rel_fitness = []
        chromosome_selected = []
        for ind in population:
            global_fitness += ind.fitness
        for ind in population:
            rel_fitness.append(float(ind.fitness)/global_fitness)
        probas = [sum(rel_fitness[:i+1]) for i in range(len(rel_fitness))]
        while len(chromosome_selected) < 2:
            random_number = random.random()
            for (i, ind) in enumerate(population):
                if random_number <= probas[i]:
                    if ind not in chromosome_selected:
                        chromosome_selected.append(ind)
                        break
        return chromosome_selected
