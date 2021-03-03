import sys
from enum import Enum

import numpy as np
from frozendict import frozendict
from tqdm import tqdm

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


def select_random_motif(s, k):
    n = len(s)
    i = np.random.randint(0, n - k + 1)
    return s[i : i + k]


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


def score_m(motifs):
    return score(motifs_to_profile(motifs))


def make_pseudo_profile(profile):
    return profile + 1


def motifs_to_profile(motifs):
    k = len(motifs[0])
    profile = np.zeros((k, len(Symbols)), dtype=int)
    for motif in motifs:
        add_to_profile(motif, profile)
    return profile


def randomized_motif_search_iter(dna, k):
    best_motifs = [select_random_motif(s, k) for s in dna]
    best_profile = motifs_to_profile(best_motifs)
    best_score = score(best_profile)

    while True:
        cur_motifs = [
            profile_most_probable(s, k, make_pseudo_profile(best_profile))
            for s in dna
        ]
        cur_profile = motifs_to_profile(cur_motifs)
        cur_score = score(cur_profile)

        if cur_score < best_score:
            best_score = cur_score
            best_profile = cur_profile
            best_motifs = cur_motifs
        else:
            return best_motifs, best_score


def randomized_motif_search(dna, k, n_iters=1000):
    best_motifs, best_score = randomized_motif_search_iter(dna, k)

    for _ in tqdm(range(n_iters)):
        cur_motifs, cur_score = randomized_motif_search_iter(dna, k)
        if cur_score < best_score:
            best_score = cur_score
            best_motifs = cur_motifs

    return best_motifs


def main():
    k, n = map(int, input().split())
    dna = []
    for _ in range(n):
        dna.append(input().strip())
    motifs = randomized_motif_search(dna, k)
    print(*motifs, sep="\n")
    # print(score_m(motifs))


if __name__ == "__main__":
    main()
