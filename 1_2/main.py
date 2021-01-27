import random

import Bio.Align.substitution_matrices as substitution_matrices
import numpy as np
from Bio import SeqIO
from Bio.pairwise2 import align, affine_penalty

_BLANK = "-"


def needleman_wunsch(s, t, cost_matrix, gap_cost=-1):
    n, m = len(s), len(t)

    res = np.zeros((n + 1, m + 1))
    res[..., 0] = np.arange(n + 1) * gap_cost
    res[0, ...] = np.arange(m + 1) * gap_cost

    for i in range(n):
        for j in range(m):
            res[i + 1, j + 1] = max(
                res[i][j] + cost_matrix[(s[i], t[j])],
                res[i + 1][j] + gap_cost,
                res[i][j + 1] + gap_cost
            )
    s_ans, t_ans = [], []
    i, j = n - 1, m - 1
    while i >= 0 and j >= 0:
        if res[i + 1, j + 1] == res[i, j + 1] + gap_cost:
            s_ans.append(s[i])
            t_ans.append(_BLANK)
            i -= 1
        if res[i + 1, j + 1] == res[i + 1, j] + gap_cost:
            t_ans.append(t[j])
            s_ans.append(_BLANK)
            j -= 1
        if res[i + 1, j + 1] == res[i, j] + cost_matrix[(s[i], t[j])]:
            s_ans.append(s[i])
            t_ans.append(t[j])
            i -= 1
            j -= 1
    while i >= 0:
        s_ans.append(s[i])
        t_ans.append(_BLANK)
        i -= 1
    while j >= 0:
        t_ans.append(t[j])
        s_ans.append(_BLANK)
        j -= 1

    return res[n][m], "".join(reversed(s_ans)), "".join(reversed(t_ans))


def affine_gap_alignment(s, t, cost_matrix, alpha, beta):
    INF = 1e9
    n, m = len(s), len(t)

    resM = np.zeros((n + 1, m + 1))
    resA = resM.copy()
    resB = resM.copy()

    resA[1:, 0] = -INF
    resA[0, 1:] = np.arange(1, m + 1) * beta + alpha

    resB[1:, 0] = np.arange(1, n + 1) * beta + alpha
    resB[0, 1:] = -INF

    resM[1:, 0] = -INF
    resM[0, 1:] = -INF

    for i in range(n):
        for j in range(m):
            resA[i + 1, j + 1] = max(
                resM[i + 1, j] + alpha + beta,
                resA[i + 1, j] + beta,
                resB[i + 1, j] + alpha + beta
            )
            resB[i + 1, j + 1] = max(
                resM[i, j + 1] + alpha + beta,
                resA[i, j + 1] + alpha + beta,
                resB[i, j + 1] + beta
            )
            resM[i + 1, j + 1] = cost_matrix[(s[i], t[j])] + max(
                resM[i, j], resA[i, j], resB[i, j]
            )

    ans = max(resM[n, m], resA[n, m], resB[n, m])
    if ans == resM[n, m]:
        cur = resM
    elif ans == resA[n, m]:
        cur = resA
    else:
        cur = resB

    s_ans, t_ans = [], []
    i, j = n - 1, m - 1
    while i >= 0 and j >= 0:
        if cur is resA:
            s_ans.append(_BLANK)
            t_ans.append(t[j])
            for mat in (resB, resM):
                if cur[i + 1, j + 1] == mat[i + 1, j] + alpha + beta:
                    cur = mat
            j -= 1
        elif cur is resB:
            s_ans.append(s[i])
            t_ans.append(_BLANK)
            for mat in (resA, resM):
                if cur[i + 1, j + 1] == mat[i, j + 1] + alpha + beta:
                    cur = mat
            i -= 1
        else:
            s_ans.append(s[i])
            t_ans.append(t[j])
            mx = max(resA[i, j], resB[i, j], resM[i, j])
            if mx == resA[i, j]:
                cur = resA
            elif mx == resB[i, j]:
                cur = resB
            i -= 1
            j -= 1

    while i >= 0:
        s_ans.append(s[i])
        t_ans.append(_BLANK)
        i -= 1
    while j >= 0:
        t_ans.append(t[j])
        s_ans.append(_BLANK)
        j -= 1

    return ans, "".join(reversed(s_ans)), "".join(reversed(t_ans))


def alignment_score(a, b, cost_matrix, gap_cost=-1):
    assert len(a) == len(b)
    score = 0
    for x, y in zip(a, b):
        if x == _BLANK or y == _BLANK:
            score += gap_cost
        else:
            score += cost_matrix[(x, y)]
    return score


def affine_alignment_score(a, b, cost_matrix, alpha, beta):
    assert len(a) == len(b)
    score = 0
    gap_x = gap_y = False
    for x, y in zip(a, b):
        if x == _BLANK:
            if not gap_x:
                gap_x = True
                gap_y = False
                score += alpha + beta
            else:
                score += beta
        elif y == _BLANK:
            if not gap_y:
                gap_y = True
                gap_x = False
                score += alpha + beta
            else:
                score += beta
        else:
            score += cost_matrix[(x, y)]
            gap_x = gap_y = False
    return score


def test_needleman_wunsch(s, t, cost_matrix, gap_cost=-1):
    score, *alignments = needleman_wunsch(s, t, cost_matrix, gap_cost)
    true_score = align.globalds(s, t, cost_matrix, gap_cost, gap_cost)[0].score
    assert score == true_score
    assert true_score == alignment_score(*alignments, cost_matrix, gap_cost)


def test_affine_gap(s, t, cost_matrix, alpha, beta):
    score, *alignments = affine_gap_alignment(s, t, cost_matrix, alpha, beta)
    test_score = affine_alignment_score(*alignments, cost_matrix, alpha, beta)
    assert score == test_score
    gap_fn = affine_penalty(alpha + beta, beta, False)
    true_score = align.globaldc(s, t, cost_matrix, gap_fn, gap_fn)[0].score
    assert score == true_score


def random_dna(length):
    return "".join(random.choices(("A", "G", "T", "C"), k=length))


def main():
    alphas = [-1, -2]
    betas = [-1, -2]
    gap_costs = [-1, -2]

    cost_matrix = substitution_matrices.load("BLOSUM62")

    # seq tests
    for fname in ("../1_1/data/gattaca.fasta", "../1_1/data/GATTACA2.fasta"):
        seqs = []
        for record in SeqIO.parse(fname, "fasta"):
            seqs.append(record.seq)
        print(fname, end=" ")
        for alpha in alphas:
            for beta in betas:
                test_affine_gap(*seqs, cost_matrix, alpha, beta)
        for gap_cost in gap_costs:
            test_needleman_wunsch(*seqs, cost_matrix, gap_cost)
        print("passed")

    # rand tests
    n_tests = 100
    max_len = 20

    for i in range(n_tests):
        s = random_dna(random.randrange(max_len - 1) + 1)
        t = random_dna(random.randrange(max_len - 1) + 1)
        for alpha in alphas:
            for beta in betas:
                test_affine_gap(s, t, cost_matrix, alpha, beta)
        for gap_cost in gap_costs:
            test_needleman_wunsch(s, t, cost_matrix, gap_cost)
    print("random tests passed")


if __name__ == "__main__":
    main()
