import sys
from enum import Enum

import numpy as np
from frozendict import frozendict
from numpy.core.numeric import Infinity

input = sys.stdin.readline


class Symbols(Enum):

    A = "A"
    C = "C"
    G = "G"
    T = "T"


class SymbolsMapping:
    idx_to_char = tuple(c.value for c in Symbols)
    char_to_idx = frozendict({k: i for i, k in enumerate(idx_to_char)})


def window(s, k):
    n = len(s)
    for i in range(n - k + 1):
        yield s[i : i + k]


def add_to_profile(s, profile):
    for i, c in enumerate(s):
        profile[i, SymbolsMapping.char_to_idx[c]] += 1
    return profile


def profile_most_probable(s, k, profile):
    best_prob = 0
    best_motif = s[:k]

    for motif in window(s, k):
        prob = 1
        for i, c in enumerate(motif):
            prob *= profile[i, SymbolsMapping.char_to_idx[c]]
        if prob > best_prob:
            best_prob = prob
            best_motif = motif
    return best_motif


def score(profile):
    return (profile.sum(axis=1) - profile.max(axis=1)).sum()


def greedy_motif_search(dna, k):
    best_score = Infinity
    best_motifs = []
    for start_motif in window(dna[0], k):
        profile = np.zeros((k, len(Symbols)), dtype=int)
        add_to_profile(start_motif, profile)
        motifs = [start_motif]
        for i in range(1, len(dna)):
            motifs.append(profile_most_probable(dna[i], k, profile))
            add_to_profile(motifs[-1], profile)
        cur_score = score(profile)
        if cur_score < best_score:
            best_score = cur_score
            best_motifs = motifs
    return best_motifs


def main():
    k, n = map(int, input().split())

    dna = []
    for _ in range(n):
        dna.append(input().strip())

    motifs = greedy_motif_search(dna, k)
    print(*motifs, sep="\n")


if __name__ == "__main__":
    main()
