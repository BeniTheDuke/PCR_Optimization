import numpy as np
from numpy.random import Generator

from CostCalculator import CostCalculator
from OptimizationResult import OptimizationResult
from Optimizer import Optimizer
from Sequence import Sequence
from Nucleotide import Nucleotide

NUM_NUCLEOTIDES = 4


class AntColonyOptimizer(Optimizer):
    """Ant System: construct full primer pairs from pheromone; deposit Q/cost on edges."""

    N_ANTS = 10
    ALPHA = 1.0
    RHO = 0.1
    Q = 1.0
    EPS = 1e-9

    @staticmethod
    def optimize(
        initial_sequences: tuple[Sequence, Sequence], max_steps: int, rng: Generator
    ) -> OptimizationResult:
        initial_cost = CostCalculator.calculate_total_cost(initial_sequences)
        result = OptimizationResult(
            (initial_sequences[0].copy(), initial_sequences[1].copy()),
            initial_cost,
            [],
        )

        tau = np.ones((2, Sequence.LENGTH, NUM_NUCLEOTIDES), dtype=np.float64)
        n_decisions = 2 * Sequence.LENGTH

        for _ in range(max_steps):
            ant_solutions: list[tuple[float, list[list[Nucleotide]]]] = []

            for _ in range(AntColonyOptimizer.N_ANTS):
                primers: list[list[Nucleotide]] = [
                    [Nucleotide.T] * Sequence.LENGTH,
                    [Nucleotide.T] * Sequence.LENGTH,
                ]

                for step in range(n_decisions):
                    if step < Sequence.LENGTH:
                        primer_idx = 0
                        pos = step
                    else:
                        primer_idx = 1
                        pos = step - Sequence.LENGTH

                    probs = tau[primer_idx, pos, :] ** AntColonyOptimizer.ALPHA
                    s = probs.sum()
                    if s <= 0 or not np.isfinite(s):
                        probs = np.ones(NUM_NUCLEOTIDES) / NUM_NUCLEOTIDES
                    else:
                        probs = probs / s

                    choice = int(rng.choice(NUM_NUCLEOTIDES, p=probs))
                    primers[primer_idx][pos] = Nucleotide(choice)

                seq0 = Sequence(primers[0])
                seq1 = Sequence(primers[1])
                cost = CostCalculator.calculate_total_cost((seq0, seq1))
                ant_solutions.append((cost, primers))

                if cost < result.best_cost:
                    result.best_cost = cost
                    result.best_sequences = (seq0.copy(), seq1.copy())

            tau *= 1.0 - AntColonyOptimizer.RHO
            for cost, primers in ant_solutions:
                deposit = AntColonyOptimizer.Q / (cost + AntColonyOptimizer.EPS)
                for step in range(n_decisions):
                    if step < Sequence.LENGTH:
                        primer_idx = 0
                        pos = step
                    else:
                        primer_idx = 1
                        pos = step - Sequence.LENGTH
                    nid = primers[primer_idx][pos].value
                    tau[primer_idx, pos, nid] += deposit

            result.cost_history.append(result.best_cost)

        return result
