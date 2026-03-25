import numpy as np
from CostCalculator import CostCalculator
from Sequence import Sequence

SEED = 42

rng = np.random.default_rng(SEED)

s1 = Sequence(Sequence._random_sequence(rng))

print(s1)
print(s1.calculate_GC_content())
print(s1.calculate_melting_temperature())
print(s1.calculate_longest_repeat())


target_melting_temperature = 50
target_gc_content = 0.5
print(CostCalculator.calculate_total_weighted_cost(s1, target_melting_temperature, target_gc_content))
