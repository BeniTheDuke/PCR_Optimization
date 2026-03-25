from numpy.random import Generator
from Nucleotide import Nucleotide


class Sequence:
    LENGTH = 20

    def __init__(self, nucleotides: list[Nucleotide]):
        if len(nucleotides) != Sequence.LENGTH:
            raise ValueError(f"Sequence not of length {Sequence.LENGTH}")

        self.nucleotides = nucleotides

    @staticmethod
    def _random_sequence(rng: Generator) -> list[Nucleotide]:
        nucleotides: list[Nucleotide] = []
        for _ in range(Sequence.LENGTH):
            value = rng.integers(0, 4)
            nucleotides.append(Nucleotide(value))

        return nucleotides

    def __str__(self):
        return "".join(n.name for n in self.nucleotides)

    def get_nucleotides_count(self) -> dict[Nucleotide, int]:
        dict = {
            Nucleotide.C: 0,
            Nucleotide.T: 0,
            Nucleotide.A: 0,
            Nucleotide.G: 0,
        }

        for n in self.nucleotides:
            dict[n] += 1

        return dict

    def calculate_GC_content(self) -> float:
        count = self.get_nucleotides_count()

        return (count[Nucleotide.G] + count[Nucleotide.C]) / Sequence.LENGTH

    def calculate_melting_temperature(self) -> float:
        count = self.get_nucleotides_count()

        return 2 * (count[Nucleotide.A] + count[Nucleotide.T]) + 4 * (
            count[Nucleotide.G] + count[Nucleotide.C]
        )

    def calculate_longest_repeat(self) -> int:
        longest_repeat = 1
        current_repeat = 1

        for i in range(1, Sequence.LENGTH):
            if self.nucleotides[i] is self.nucleotides[i - 1]:
                current_repeat += 1
            else:
                current_repeat = 1

            if current_repeat > longest_repeat:
                longest_repeat = current_repeat

        return longest_repeat

    def calculate_hairpin_score(self, min_stem: int = 3, min_loop: int = 3) -> int:
        best_score = 0

        for stem_len in range(min_stem, Sequence.LENGTH // 2 + 1):
            for i in range(Sequence.LENGTH - stem_len + 1):
                left_stem = self.nucleotides[i : i + stem_len]

                for j in range(i + stem_len + min_loop, Sequence.LENGTH - stem_len + 1):
                    right_stem = self.nucleotides[j : j + stem_len]

                    if Sequence.is_reverse_complement(left_stem, right_stem):
                        best_score = max(best_score, stem_len)

        return best_score

    @staticmethod
    def is_reverse_complement(seq_1: list[Nucleotide], seq_2: list[Nucleotide]):
        if len(seq_1) != len(seq_2):
            return False

        seq_length = len(seq_1)

        for i in range(seq_length):
            if not seq_1[i].complement != seq_2[seq_length - 1 - i]:
                return False

        return True
