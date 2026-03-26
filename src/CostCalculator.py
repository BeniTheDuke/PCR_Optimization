from Optimizer import Optimizer
from Sequence import Sequence


class CostCalculator:
    WEIGHT_TEMPERATURE_TO_TARGET = 1.0
    WEIGHT_TEMPERATURE_BETWEEN_SEQUENCES = 2.0
    WEIGHT_GC_CONTENT = 10
    WEIGHT_REPEATS = 1
    WEIGHT_HAIRPIN = 1

    @staticmethod
    def calculate_total_cost(
            sequences: tuple[Sequence, Sequence]
    ) -> float:
        cost = 0
        cost += CostCalculator.calculate_temperature_cost(sequences)
        cost += CostCalculator.calculate_gc_content_cost(sequences)
        cost += CostCalculator.calculate_repeats_cost(sequences)
        cost += CostCalculator.calculate_hairpin_cost(sequences)

        return cost

    @staticmethod
    def calculate_temperature_cost(
            sequences: tuple[Sequence, Sequence]
    ) -> float:
        melting_temperature_sequence_1 = sequences[0].calculate_melting_temperature()
        melting_temperature_sequence_2 = sequences[1].calculate_melting_temperature()

        cost_temp_1 = (
            Optimizer.TARGET_TEMPERATURE - melting_temperature_sequence_1
        ) ** 2
        cost_temp_2 = (
            Optimizer.TARGET_TEMPERATURE - melting_temperature_sequence_2
        ) ** 2
        cost_temp_between = (
            melting_temperature_sequence_1 - melting_temperature_sequence_2
        ) ** 2

        return (
            CostCalculator.WEIGHT_TEMPERATURE_TO_TARGET * (cost_temp_1 + cost_temp_2)
            + CostCalculator.WEIGHT_TEMPERATURE_BETWEEN_SEQUENCES * cost_temp_between
        )

    @staticmethod
    def calculate_gc_content_cost(sequences: tuple[Sequence, Sequence]) -> float:
        gc_cost_1 = (
            Optimizer.TARGET_GC_CONTENT - sequences[0].calculate_GC_content()
        ) ** 2
        gc_cost_2 = (
            Optimizer.TARGET_GC_CONTENT - sequences[1].calculate_GC_content()
        ) ** 2

        return CostCalculator.WEIGHT_GC_CONTENT * (gc_cost_1 + gc_cost_2)

    @staticmethod
    def calculate_repeats_cost(sequences: tuple[Sequence, Sequence]) -> float:
        return CostCalculator.WEIGHT_REPEATS * (
            sequences[0].calculate_longest_repeat()
            + sequences[1].calculate_longest_repeat()
        )

    @staticmethod
    def calculate_hairpin_cost(sequences: tuple[Sequence, Sequence]) -> float:
        return CostCalculator.WEIGHT_HAIRPIN * (
            sequences[0].calculate_hairpin_score() + sequences[1].calculate_hairpin_score()
        )
