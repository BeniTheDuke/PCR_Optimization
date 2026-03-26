from numpy.random import Generator
import sys

from CostCalculator import CostCalculator
from OptimizationResult import OptimizationResult
from Optimizer import Optimizer
from Sequence import Sequence
from Nucleotide import Nucleotide


class MultipleRandomOptimizer(Optimizer):
    CHANGES_PER_STEP = 6

    @staticmethod
    def optimize(
        initial_sequences: tuple[Sequence, Sequence], max_steps: int, rng: Generator
    ) -> OptimizationResult:

        initial_score = CostCalculator.calculate_total_cost(initial_sequences)
        optimizationResult = OptimizationResult(initial_sequences, initial_score, [])

        for _ in range(max_steps):
            new_sequences = (
                optimizationResult.best_sequences[0].copy(),
                optimizationResult.best_sequences[1].copy(),
            )

            for _ in range(MultipleRandomOptimizer.CHANGES_PER_STEP):
                sequence_to_change = rng.integers( 0, 2)  # which one do we want to change?
                change_position = rng.integers(0, Sequence.LENGTH)
                new_nucleotide = Nucleotide(rng.integers(0, len(Nucleotide)))
                new_sequences[sequence_to_change].nucleotides[change_position] = (
                    new_nucleotide
                )

            new_cost = CostCalculator.calculate_total_cost(new_sequences)

            if new_cost <= optimizationResult.best_cost:
                optimizationResult.best_cost = new_cost
                optimizationResult.best_sequences = new_sequences

            optimizationResult.cost_history.append(new_cost)

        return optimizationResult
