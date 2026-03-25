import numpy as np
import hairpin_dimer_scores
from Sequence import Sequence
from CostCalculator import CostCalculator

SEED = 42

target_melting_temperature = 50
target_gc_content = 0.5


rng = np.random.default_rng(SEED)

sequences: list[Sequence] = []
for i in range(3):
    sequences.append(Sequence(Sequence._random_sequence(rng)))

for sequence in sequences:
    print("----------------")
    print(sequence)
    print(f"GC Content: {sequence.calculate_GC_content()}")
    print(f"Melting temperature: {sequence.calculate_melting_temperature()}")
    print(f"Longest repeat: {sequence.calculate_longest_repeat()}")
    print(f"Hairspin score: {sequence.calculate_hairpin_score()}")
    cost = CostCalculator.calculate_total_weighted_cost(sequence, target_melting_temperature, target_gc_content)
    print(f"Cost: {cost}")
    print("----------------")


