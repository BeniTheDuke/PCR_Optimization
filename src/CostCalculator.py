from Sequence import Sequence


class CostCalculator:
    WEIGHT_TEMPERATURE = 1
    WEIGHT_GC_CONTENT = 10
    WEIGHT_REPEATS = 1

    @staticmethod
    def calculate_total_weighted_cost(
        sequence: Sequence, target_melting_temperature, target_gc_content
    ) -> float:
        cost_temperature = (
            CostCalculator.WEIGHT_TEMPERATURE
            * CostCalculator.calculate_temperature_cost(sequence, target_melting_temperature)
        )
        cost_gc_content = (
            CostCalculator.WEIGHT_GC_CONTENT
            * CostCalculator.calculate_gc_content_cost(sequence, target_gc_content)
        )
        cost_repeats = (
            CostCalculator.WEIGHT_REPEATS * sequence.calculate_longest_repeat()
        )

        return cost_temperature + cost_gc_content + cost_repeats

    @staticmethod
    def calculate_temperature_cost(sequence: Sequence, target_melting_temperature) -> float:
        return(target_melting_temperature - sequence.calculate_melting_temperature()) ** 2 

    @staticmethod
    def calculate_gc_content_cost(sequence: Sequence, target_gc_content) -> float:
        return (target_gc_content - sequence.calculate_GC_content()) ** 2

