from abc import ABC, abstractmethod
from numpy.random import Generator

from Sequence import Sequence
from OptimizationResult import OptimizationResult


class Optimizer(ABC):
    TARGET_TEMPERATURE = 60
    TARGET_GC_CONTENT = 0.5

    @staticmethod
    @abstractmethod
    def optimize(
        initial_sequences: tuple[Sequence, Sequence],
        max_steps: int,
        rng: Generator,
    ) -> OptimizationResult:
        pass
