from __future__ import annotations
from enum import Enum


class Nucleotide(Enum):
    T = 0
    C = 1
    A = 2
    G = 3

    def __str__(self):
        return self.name

    def complement(self) -> Nucleotide:
        complements = {
            Nucleotide.A: Nucleotide.T,
            Nucleotide.T: Nucleotide.A,
            Nucleotide.C: Nucleotide.G,
            Nucleotide.G: Nucleotide.C,
        }

        return complements[self]
