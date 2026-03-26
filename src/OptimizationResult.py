from dataclasses import dataclass
from Sequence import Sequence

@dataclass
class OptimizationResult:
    best_sequences: tuple[Sequence, Sequence] # forward, backward
    best_cost: float
    cost_history: list[float] 
