from numpy.random import Generator
import numpy as np
import matplotlib.pyplot as plt

from CostCalculator import CostCalculator
from OptimizationResult import OptimizationResult
from Optimizer import Optimizer
from Sequence import Sequence
from Nucleotide import Nucleotide


class GeneticAlgorithmOptimizer(Optimizer):
    @staticmethod
    def random_primer_pair(
        rng:Generator):
        sequence_1 = Sequence(Sequence._random_sequence(rng))
        sequence_2 = Sequence(Sequence._random_sequence(rng))

        return sequence_1, sequence_2

    @staticmethod
    def initialize_population(
            population_size:int,
            rng:Generator):
        
        return [GeneticAlgorithmOptimizer.random_primer_pair(rng) 
                for _ in range(population_size)]


    @staticmethod
    def fitness(
            primer_pair):
        
        return CostCalculator.calculate_total_cost(primer_pair)


    @staticmethod
    def tournament_selection(
        population, fitnesses, rng:Generator, k=3):
        selected = []
        for _ in range(len(population)):
            idx = rng.choice(len(population), k, replace=False)
            best_idx = idx[np.argmin([fitnesses[i] for i in idx])]
            selected.append(population[best_idx])

        return selected
    

    @staticmethod
    def crossover_primer_pair(
        parent_primer_pair1, parent_primer_pair2,rng:Generator):
        cut = rng.integers(1, Sequence.LENGTH)
        
        #create offspring from parents primer
        child1_seq1 = parent_primer_pair1[0].copy()
        child1_seq2 = parent_primer_pair1[1].copy()
        child2_seq1 = parent_primer_pair2[0].copy()
        child2_seq2 = parent_primer_pair2[1].copy()

        # Swap forward primer
        child1_seq1.nucleotides[:cut] = parent_primer_pair2[0].nucleotides[:cut]
        child2_seq1.nucleotides[:cut] = parent_primer_pair1[0].nucleotides[:cut]

        # Swap reverse primer
        child1_seq2.nucleotides[:cut] = parent_primer_pair2[1].nucleotides[:cut]
        child2_seq2.nucleotides[:cut] = parent_primer_pair1[1].nucleotides[:cut]

        return (child1_seq1, child1_seq2), (child2_seq1, child2_seq2)
    

    @staticmethod
    def mutate_primer_pair(
        primer_pair, rng:Generator, mutation_rate = 1/(Sequence.LENGTH*2)): # 1/L Heuristic (Gabriela Ochoa, 2002)

        sequence_1 = primer_pair[0].copy() 
        sequence_2 = primer_pair[1].copy()

        for sequence in (sequence_1, sequence_2):
            for i in range(Sequence.LENGTH):
                if rng.random() < mutation_rate:
                    sequence.nucleotides[i] = Nucleotide(rng.integers(0, len(Nucleotide)))
        
        return (sequence_1, sequence_2)
    

    @staticmethod
    def optimize(
        initial_sequences, max_steps: int, rng:Generator, 
        population_size=4,
        k=3
        )-> OptimizationResult:

        population = [initial_sequences] + GeneticAlgorithmOptimizer.initialize_population(
            population_size - 1, rng)

        initial_score = CostCalculator.calculate_total_cost(initial_sequences)
        result = OptimizationResult(initial_sequences, initial_score, [])

        best_cost = initial_score
        best_primer_pair = initial_sequences

        generations = max_steps

        for gen in range(generations):

            fitnesses = [
                GeneticAlgorithmOptimizer.fitness(primer_pair) for primer_pair in population
            ]

            best_idx = np.argmin(fitnesses)
            gen_best_cost = fitnesses[best_idx]

            if gen_best_cost < best_cost:
                best_cost = gen_best_cost
                best_primer_pair = (
                    population[best_idx][0].copy(),
                    population[best_idx][1].copy()
                )

            # tracking
            result.cost_history.append(best_cost)
            result.best_cost = best_cost
            result.best_sequences = best_primer_pair

            # selection
            parents = GeneticAlgorithmOptimizer.tournament_selection(
                population, fitnesses, rng, k
            )

            rng.shuffle(parents)

            next_population = []

            # crossover + mutation
            for i in range(0, len(parents) - 1, 2):
                p1 = parents[i]
                p2 = parents[i + 1]

                child1, child2 = GeneticAlgorithmOptimizer.crossover_primer_pair(p1, p2, rng)

                next_population.append(
                    GeneticAlgorithmOptimizer.mutate_primer_pair(child1, rng)
                )
                next_population.append(
                    GeneticAlgorithmOptimizer.mutate_primer_pair(child2, rng)
                )

            # ensure correct population size
            if len(next_population) > population_size:
                next_population = next_population[:population_size]

            while len(next_population) < population_size:
                next_population.append(
                    GeneticAlgorithmOptimizer.random_primer_pair(rng)
                )

            # elitism
            fitnesses_next = [
                GeneticAlgorithmOptimizer.fitness(pair)
                for pair in next_population
            ]

            worst_idx = np.argmax(fitnesses_next)

            next_population[worst_idx] = (
                best_primer_pair[0].copy(),
                best_primer_pair[1].copy()
            )
    

        
            population = next_population

        return result
  