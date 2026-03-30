import numpy as np
from numpy.random import default_rng

from CostCalculator import CostCalculator
from Sequence import Sequence
from Nucleotide import Nucleotide

NUM_NUCLEOTIDES = 4

def ant_colony_optimize(initial_sequences, max_steps, n_ants):
    """
        Perform Ant Colony Optimization (ACO) to minimize the cost of two sequences.

        The algorithm simulates a colony of ants constructing candidate sequences
        guided by pheromone trails. Pheromones are updated based on the quality
        (cost) of solutions, allowing the colony to converge toward optimal sequences.

        Parameters
        ----------
        initial_sequences : tuple of Sequence
            Tuple containing two initial sequences to start optimization.
        max_steps : int
            Maximum number of iterations for the optimization process.
        n_ants : int
            Number of ants (candidate solutions) generated per iteration.

        Returns
        -------
        best_sequences : tuple of Sequence
            The best sequences found during the optimization.
        best_cost : float
            Cost associated with the best sequences.
        cost_history : list of float
            History of best costs recorded at each iteration.
        """

    rng = default_rng()


    best_sequences = (initial_sequences[0].copy(), initial_sequences[1].copy())
    best_cost = CostCalculator.calculate_total_cost(best_sequences)
    cost_history = []

    tau = np.ones((2, Sequence.LENGTH, NUM_NUCLEOTIDES))

    ALPHA = 1.0
    RHO = 0.1
    Q = 1.0
    EPS = 1e-9

    for step in range(max_steps):
        all_solutions = []

        for ant in range(n_ants):
            primers = [
                [Nucleotide.T] * Sequence.LENGTH,
                [Nucleotide.T] * Sequence.LENGTH
            ]

            for i in range(2 * Sequence.LENGTH):
                seq_idx = 0 if i < Sequence.LENGTH else 1
                pos = i if i < Sequence.LENGTH else i - Sequence.LENGTH

                probs = tau[seq_idx][pos] ** ALPHA
                total = probs.sum()
                if total <= 0 or not np.isfinite(total):
                    probs = np.ones(NUM_NUCLEOTIDES) / NUM_NUCLEOTIDES
                else:
                    probs = probs / total

                choice = rng.choice(NUM_NUCLEOTIDES, p=probs)
                primers[seq_idx][pos] = Nucleotide(int(choice))

            seq0 = Sequence(primers[0])
            seq1 = Sequence(primers[1])
            cost = CostCalculator.calculate_total_cost((seq0, seq1))
            all_solutions.append((cost, (seq0, seq1)))

            #if cost < best_cost:            #can be unchecked or not just helps with visualisation if unchecked
            best_cost = cost
            best_sequences = (seq0.copy(), seq1.copy())

        tau *= (1 - RHO)

        for cost, seqs in all_solutions:
            deposit = Q / (cost + EPS)
            for seq_idx in range(2):
                for pos in range(Sequence.LENGTH):
                    nid = seqs[seq_idx].nucleotides[pos].value
                    tau[seq_idx][pos][nid] += deposit

        cost_history.append(best_cost)

    return best_sequences, best_cost, cost_history


# Test
if __name__ == "__main__":
    rng = default_rng(42)

    # Initial sequences: all T
    seq_a = Sequence([Nucleotide.T] * Sequence.LENGTH)
    seq_b = Sequence([Nucleotide.T] * Sequence.LENGTH)

    best_sequences, best_cost, history = ant_colony_optimize((seq_a, seq_b), max_steps=100, n_ants=10)

    print("Best cost:", best_cost)
    print("Cost history:", history)
    print("Best sequences:", best_sequences)
