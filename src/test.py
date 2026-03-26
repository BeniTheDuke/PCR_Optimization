import numpy as np
from MultipleRandomOptimizer import MultipleRandomOptimizer
from Optimizer import Optimizer
from RandomOptimizer import RandomOptimizer
from Sequence import Sequence
from CostCalculator import CostCalculator

SEED = 42


rng = np.random.default_rng(SEED)

sequence_1 = Sequence(Sequence._random_sequence(rng))
sequence_2 = Sequence(Sequence._random_sequence(rng))


print("Initial Sequences")
print(sequence_1)
print(sequence_2)

print("Optimization:")
result = RandomOptimizer.optimize((sequence_1, sequence_2), 1000, rng)

print(result.cost_history)
print(result.best_sequences[0])
print(result.best_sequences[1])
