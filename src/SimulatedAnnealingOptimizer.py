import math
from numpy.random import Generator

from CostCalculator import CostCalculator
from OptimizationResult import OptimizationResult
from Optimizer import Optimizer
from Sequence import Sequence
from Nucleotide import Nucleotide


class SimulatedAnnealingOptimizer(Optimizer):
    INITIAL_TEMPERATURE = 100.0
    COOLING_RATE = 0.995
    MIN_TEMPERATURE = 0.01

    @staticmethod
    def optimize(
        initial_sequences: tuple[Sequence, Sequence], max_steps: int, rng: Generator
    ) -> OptimizationResult:

        current_sequences = (initial_sequences[0].copy(), initial_sequences[1].copy())
        current_cost = CostCalculator.calculate_total_cost(current_sequences)

        best_sequences = (current_sequences[0].copy(), current_sequences[1].copy())
        best_cost = current_cost
        cost_history = []

        temperature = SimulatedAnnealingOptimizer.INITIAL_TEMPERATURE

        for _ in range(max_steps):
            # copy both primers
            new_sequences = (
                current_sequences[0].copy(),
                current_sequences[1].copy(),
            )

            # mutate one random position in one random primer
            sequence_to_change = rng.integers(0, 2)
            change_position = rng.integers(0, Sequence.LENGTH)
            new_nucleotide = Nucleotide(rng.integers(0, len(Nucleotide)))
            new_sequences[sequence_to_change].nucleotides[change_position] = new_nucleotide

            new_cost = CostCalculator.calculate_total_cost(new_sequences)
            delta = new_cost - current_cost

            # accept if better, or probabilistically if worse
            if delta <= 0 or rng.random() < math.exp(-delta / temperature):
                current_sequences = new_sequences
                current_cost = new_cost

            # track global best
            if current_cost < best_cost:
                best_cost = current_cost
                best_sequences = (current_sequences[0].copy(), current_sequences[1].copy())

            cost_history.append(best_cost)

            # cool down
            temperature = max(
                SimulatedAnnealingOptimizer.MIN_TEMPERATURE,
                temperature * SimulatedAnnealingOptimizer.COOLING_RATE
            )

        return OptimizationResult(best_sequences, best_cost, cost_history)