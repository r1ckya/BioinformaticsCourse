from enum import Enum

import numpy as np
from frozendict import frozendict


class Symbols(Enum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"
    Gap = "-"


class SymbolsMapping:
    idx_to_char = tuple(c.value for c in Symbols)
    char_to_idx = frozendict({k: i for i, k in enumerate(idx_to_char)})


def symbol_to_profile(c):
    profile = np.zeros(len(Symbols))
    profile[SymbolsMapping.char_to_idx[c]] = 1.0
    return profile


class AlignmentItem:
    def __init__(self, seqs, idxs, profile, consensus):
        self.seqs = seqs
        self.idxs = idxs
        self.profile = profile
        self.consensus = consensus
        self.seqs_length = len(self.seqs[0])

    @classmethod
    def from_seq(cls, seq, idx=0):
        seqs = [list(seq)]
        profile = [symbol_to_profile(c) for c in seq]
        consensus = list(seq)
        return cls(seqs, [idx], profile, consensus)

    def merge(self, other):
        assert self.seqs_length == other.seqs_length
        n1, n2 = len(self.seqs), len(other.seqs)
        nn = n1 + n2
        self.profile = [(p1 * n1 + p2 * n2) / nn
                        for p1, p2 in zip(self.profile, other.profile)]
        self.seqs.extend(other.seqs)
        self.idxs.extend(other.idxs)
        self.consensus = [SymbolsMapping.idx_to_char[p.argmax()] for p in
                          self.profile]

    def reorder_seqs(self):
        seqs = [0] * len(self.idxs)
        for i, seq in zip(self.idxs, self.seqs):
            seqs[i] = seq
        self.seqs = seqs
        self.idxs = [i for i in range(len(self.seqs))]

    def __repr__(self):
        return repr(self.seqs)
