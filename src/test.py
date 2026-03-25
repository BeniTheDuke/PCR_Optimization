import numpy as np
from MultipleRandomOptimizer import MultipleRandomOptimizer
from Optimizer import Optimizer
from RandomOptimizer import RandomOptimizer
from Sequence import Sequence
from CostCalculator import CostCalculator

SEED = 42

target_melting_temperature = 50
target_gc_content = 0.5


rng = np.random.default_rng(SEED)

sequence = Sequence(Sequence._random_sequence(rng))


print("Initial Sequence")
sequence.print()

print("Optimization:")
result = MultipleRandomOptimizer.optimize(sequence, 1000, rng)

print(result.cost_history)
result.best_sequence.print()

