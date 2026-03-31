import numpy as np
from MultipleRandomOptimizer import MultipleRandomOptimizer
from SimulatedAnnealingOptimizer import SimulatedAnnealingOptimizer
from GeneticAlgorithmOptimizer import GeneticAlgorithmOptimizer
from AntColonyOptimizer import AntColonyOptimizer
from Optimizer import Optimizer
from RandomOptimizer import RandomOptimizer
from Sequence import Sequence
from CostCalculator import CostCalculator
from PerformanceMeasurements import PerformanceMeasurements

SEED = 42


rng = np.random.default_rng(SEED)

sequence_1 = Sequence(Sequence.random_sequence(rng))
sequence_2 = Sequence(Sequence.random_sequence(rng))

run_random = True
# run_random = False
run_sa = True
# run_sa = False
run_ga = True
# run_ga = False
run_ant = True
# run_ant = False

print("Initial Sequences")
print("Forward Primer: ", sequence_1)
print("Reverse Primer: ",sequence_2)

### Random Optimizer
if run_random:
    print("Random Optimization:")
    rand_result, rand_performance = RandomOptimizer.optimize_with_performance_measurements((sequence_1, sequence_2), 1000, rng)
    print(rand_result.cost_history)
    print("Rand Final Forward primer: ", rand_result.best_sequences[0])
    print("Rand Final Reverse primer: ", rand_result.best_sequences[1])
    print(f"Rand Best cost: {rand_result.best_cost}")
    print(f"Rand Duration: {rand_performance.runtime_seconds}")
    print(f"Rand Peak memory usage: {rand_performance.peak_memory_bytes}")

### Simulated Annealing
if run_sa:
    print("Simulated Annealing Optimization:" )
    sa_result, sa_performance = SimulatedAnnealingOptimizer.optimize_with_performance_measurements((sequence_1, sequence_2), 1000, rng)
    print(sa_result.cost_history)
    print("SA Final Forward primer: ", sa_result.best_sequences[0])
    print("SA Final Reverse primer: ", sa_result.best_sequences[1])
    print(f"SA Best cost: {sa_result.best_cost}")
    print(f"SA Duration: {sa_performance.runtime_seconds}")
    print(f"SA Peak memory usage: {sa_performance.peak_memory_bytes}")

### Genetic Algorithm
if run_ga:
    print("Genetic Algorithm Optimization:" )
    ga_result, ga_performance = GeneticAlgorithmOptimizer.optimize_with_performance_measurements((sequence_1, sequence_2), 1000, rng)
    print(ga_result.cost_history)
    print("GA Final Forward primer: ", ga_result.best_sequences[0])
    print("GA Final Reverse primer: ", ga_result.best_sequences[1])
    print(f"GA Best cost: {ga_result.best_cost}")
    print(f"GA Duration: {ga_performance.runtime_seconds}")
    print(f"GA Peak memory usage: {ga_performance.peak_memory_bytes}")

### Ant Swarm
if run_ant:
    print("Ant Swarm Optimization:" )
    ant_result, ant_performance = AntColonyOptimizer.optimize_with_performance_measurements((sequence_1, sequence_2), 1000, rng)
    print(ant_result.cost_history)
    print("Ant Final Forward primer: ", ant_result.best_sequences[0])
    print("Ant Final Reverse primer: ", ant_result.best_sequences[1])
    print(f"Ant Best cost: {ant_result.best_cost}")
    print(f"Ant Duration: {ant_performance.runtime_seconds}")
    print(f"Ant Peak memory usage: {ant_performance.peak_memory_bytes}")