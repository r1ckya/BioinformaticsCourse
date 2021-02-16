import numpy as np

from alignment_item import Symbols, symbol_to_profile, AlignmentItem


def scorer(del_cost, ins_cost, match_cost, mismatch_cost):
    def inner(a, b):
        if a == Symbols.Gap.value:
            return ins_cost
        elif b == Symbols.Gap.value:
            return del_cost
        elif a == b:
            return match_cost
        else:
            return mismatch_cost

    return inner


def append_symbol(i, consensus, seqs, profile, align_profile):
    consensus.append(align_profile.consensus[i])
    for j in range(len(align_profile.seqs)):
        seqs[j].append(align_profile.seqs[j][i])
    profile.append(align_profile.profile[i])


def append_gap(consensus, seqs, profile):
    consensus.append(Symbols.Gap.value)
    for j in range(len(seqs)):
        seqs[j].append(Symbols.Gap.value)
    profile.append(symbol_to_profile(Symbols.Gap.value))


def reverse_alignment(consensus, seqs, profile):
    consensus.reverse()
    for i in range(len(seqs)):
        seqs[i].reverse()
    profile.reverse()


def restore_alignment(s_align_profile, t_align_profile, res,
                      del_cost, ins_cost, score_fn):
    s = s_align_profile.consensus
    t = t_align_profile.consensus
    n, m = len(s), len(t)

    s_consensus = []
    s_profile = []
    s_seqs = [[] for _ in range(len(s_align_profile.seqs))]
    t_consensus = []
    t_profile = []
    t_seqs = [[] for _ in range(len(t_align_profile.seqs))]

    i, j = n - 1, m - 1
    while i >= 0 and j >= 0:
        if res[i + 1, j + 1] == res[i, j + 1] + del_cost:
            append_symbol(i, s_consensus, s_seqs, s_profile, s_align_profile)
            append_gap(t_consensus, t_seqs, t_profile)
            i -= 1
        elif res[i + 1, j + 1] == res[i + 1, j] + ins_cost:
            append_gap(s_consensus, s_seqs, s_profile)
            append_symbol(j, t_consensus, t_seqs, t_profile, t_align_profile)
            j -= 1
        elif res[i + 1, j + 1] == res[i, j] + score_fn(s[i], t[j]):
            append_symbol(i, s_consensus, s_seqs, s_profile, s_align_profile)
            append_symbol(j, t_consensus, t_seqs, t_profile, t_align_profile)
            i -= 1
            j -= 1
    while i >= 0:
        append_symbol(i, s_consensus, s_seqs, s_profile, s_align_profile)
        append_gap(t_consensus, t_seqs, t_profile)
        i -= 1
    while j >= 0:
        append_gap(s_consensus, s_seqs, s_profile)
        append_symbol(j, t_consensus, t_seqs, t_profile, t_align_profile)
        j -= 1

    reverse_alignment(s_consensus, s_seqs, s_profile)
    reverse_alignment(t_consensus, t_seqs, t_profile)
    s_idxs = s_align_profile.idxs.copy()
    t_idxs = t_align_profile.idxs.copy()
    new_alignment = AlignmentItem(s_seqs, s_idxs, s_profile, s_consensus)
    new_alignment.merge(
        AlignmentItem(t_seqs, t_idxs, t_profile, t_consensus)
    )
    return new_alignment


def needleman_wunsch(s_align_profile, t_align_profile,
                     del_cost, ins_cost,
                     match_cost, mismatch_cost):
    s = s_align_profile.consensus
    t = t_align_profile.consensus
    score_fn = scorer(del_cost, ins_cost, match_cost, mismatch_cost)
    n, m = len(s), len(t)

    res = np.zeros((n + 1, m + 1))
    res[..., 0] = np.arange(n + 1) * del_cost
    res[0, ...] = np.arange(m + 1) * ins_cost

    for i in range(n):
        for j in range(m):
            res[i + 1, j + 1] = max(
                res[i][j] + score_fn(s[i], t[j]),
                res[i + 1][j] + ins_cost,
                res[i][j + 1] + del_cost
            )

    new_alignment = restore_alignment(s_align_profile, t_align_profile,
                                      res, del_cost, ins_cost, score_fn)
    return res[n][m], new_alignment


def needleman_wunsch_score_only(s_align_profile, t_align_profile,
                                del_cost, ins_cost,
                                match_cost, mismatch_cost):
    s = s_align_profile.consensus
    t = t_align_profile.consensus
    score_fn = scorer(del_cost, ins_cost, match_cost, mismatch_cost)
    n, m = len(s), len(t)

    if n < m:
        n, m = m, n
        s, t = t, s

    res = np.zeros((2, m + 1))
    res[0, ...] = np.arange(m + 1) * ins_cost

    k = 0
    for i in range(n):
        k ^= 1
        res[k, 0] = (i + 1) * del_cost
        for j in range(m):
            res[k, j + 1] = max(
                res[k ^ 1][j] + score_fn(s[i], t[j]),
                res[k][j] + ins_cost,
                res[k ^ 1][j + 1] + del_cost
            )

    return res[k][m]
