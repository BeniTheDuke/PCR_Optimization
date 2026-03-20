import numpy as np
from Sequence import Sequence

SEED = 42

rng = np.random.default_rng(SEED)

s1 = Sequence(Sequence._random_sequence(rng))

print(s1)
print(s1.calculate_GC_content())
print(s1.calculate_melting_temperature())
print(s1.calculate_longest_repeat())
