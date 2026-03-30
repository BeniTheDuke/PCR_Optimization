from abc import ABC, abstractmethod
from numpy.random import Generator
import tracemalloc
from time import perf_counter

from Sequence import Sequence
from OptimizationResult import OptimizationResult
from PerformanceMeasurements import PerformanceMeasurements


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


    @classmethod
    def optimize_with_performance_measurements(
        cls,
        initial_sequences: tuple[Sequence, Sequence],
        max_steps: int,
        rng: Generator,
    ) -> tuple[OptimizationResult, PerformanceMeasurements]:
            tracemalloc.start()
            start_time = perf_counter()

            optimization_result = cls.optimize(initial_sequences, max_steps, rng)

            end_time = perf_counter()
            _, peak_memory = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            performance_measurements = PerformanceMeasurements(end_time - start_time, peak_memory)

            return (optimization_result, performance_measurements)



