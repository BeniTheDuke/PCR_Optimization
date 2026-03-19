from enum import Enum


class Nucleotide(Enum):
    T = 0
    C = 1
    A = 2
    G = 3

    def __str__(self):
        return self.name
