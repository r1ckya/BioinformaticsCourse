import sys
from enum import Enum

import numpy as np


class Letters(Enum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"

    def to_tuple():
        return tuple(c.value for c in Letters)


def random_seq(vocab, length):
    return "".join(np.random.choice(vocab, size=length))


def random_dna(length):
    return random_seq(Letters.to_tuple(), length)


def random_read(dna, read_len):
    i = np.random.randint(0, len(dna) - read_len + 1)
    return dna[i : i + read_len]


def full_coverage(dna, read_len):
    for i in range(0, len(dna) - read_len + 1):
        yield dna[i : i + read_len]


def main():
    dna_len = 1000
    read_len = 150
    n_reads = 1000

    dna = random_dna(dna_len)
    with open("dna.txt", "w") as f:
        print(dna, file=f)

    for i in range(n_reads):
        read = random_read(dna, read_len)
        print(f"@{i}", file=sys.stdout)
        print(read, file=sys.stdout)
        print("+", file=sys.stdout)
        print("I" * read_len, file=sys.stdout)

        # for read in full_coverage(dna, read_len):
        #     print(f"@{i}", file=f)
        #     print(read, file=f)
        #     print("+", file=f)
        #     print("I" * read_len, file=f)


if __name__ == "__main__":
    main()
