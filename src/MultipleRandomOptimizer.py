from numpy.random import Generator
import sys

from CostCalculator import CostCalculator
from OptimizationResult import OptimizationResult
from Optimizer import Optimizer
from Sequence import Sequence
from Nucleotide import Nucleotide


class MultipleRandomOptimizer(Optimizer):
    CHANGES_PER_STEP = 10

    @staticmethod
    def optimize(
        initial_sequence: Sequence, max_steps: int, rng: Generator
    ) -> OptimizationResult:
        optimizationResult = OptimizationResult(initial_sequence, sys.maxsize, [])

        for _ in range(max_steps):
            nucleotides_copy = optimizationResult.best_sequence.nucleotides.copy()

            for _ in range(MultipleRandomOptimizer.CHANGES_PER_STEP):
                change_position = rng.integers(0, Sequence.LENGTH)
                new_nucleotide = Nucleotide(rng.integers(0, len(Nucleotide)))
                nucleotides_copy[change_position] = new_nucleotide

            new_sequence = Sequence(nucleotides_copy)
            new_cost = CostCalculator.calculate_total_weighted_cost(
                new_sequence, Optimizer.TARGET_TEMPERATURE, Optimizer.TARGET_GC_CONTENT
            )

            if new_cost <= optimizationResult.best_cost:
                optimizationResult.best_cost = new_cost
                optimizationResult.best_sequence = new_sequence

            optimizationResult.cost_history.append(optimizationResult.best_cost)

        return optimizationResult
