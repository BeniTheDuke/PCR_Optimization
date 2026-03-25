from dataclasses import dataclass
from Sequence import Sequence

@dataclass
class OptimizationResult:
    best_sequence: Sequence
    best_cost: float
    cost_history: list[float] 
