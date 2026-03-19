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
