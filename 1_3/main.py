import random

import numpy as np
from Bio.pairwise2 import align

_BLANK = "-"


def scorer(match_cost, mismatch_cost):
    return lambda a, b: match_cost if a == b else mismatch_cost


def needleman_wunsch(s, t, del_cost, ins_cost, match_cost, mismatch_cost):
    score_fn = scorer(match_cost, mismatch_cost)
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
    s_ans, t_ans = [], []
    i, j = n - 1, m - 1
    while i >= 0 and j >= 0:
        if res[i + 1, j + 1] == res[i, j + 1] + del_cost:
            s_ans.append(s[i])
            t_ans.append(_BLANK)
            i -= 1
        elif res[i + 1, j + 1] == res[i + 1, j] + ins_cost:
            t_ans.append(t[j])
            s_ans.append(_BLANK)
            j -= 1
        elif res[i + 1, j + 1] == res[i, j] + score_fn(s[i], t[j]):
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

    return res[n][m], s_ans[::-1], t_ans[::-1]


def hirschberg_score(s, t, del_cost, ins_cost, match_cost, mismatch_cost):
    n, m = len(s), len(t)
    res = np.array([
        [i * ins_cost for i in range(m + 1)],
        [0] * (m + 1)
    ])
    k = 0
    score_fn = scorer(match_cost, mismatch_cost)
    for i in range(n):
        k ^= 1
        res[k][0] = (i + 1) * del_cost
        for j in range(m):
            res[k][j + 1] = max(
                res[k ^ 1][j] + score_fn(s[i], t[j]),
                res[k][j] + ins_cost,
                res[k ^ 1][j + 1] + del_cost
            )
    return res[k]


def hirschberg(s, t, del_cost, ins_cost, match_cost, mismatch_cost,
               vis=False, depth=0):
    n, m = len(s), len(t)
    if n < 2 or m < 2:
        if vis:
            print("  " * depth, f"({s}, {t})")
        return needleman_wunsch(s, t, del_cost, ins_cost, match_cost,
                                mismatch_cost)
    else:
        mid = n // 2
        score_a = hirschberg_score(s[:mid], t, del_cost, ins_cost,
                                   match_cost, mismatch_cost)
        score_b = hirschberg_score(s[mid:][::-1], t[::-1], del_cost, ins_cost,
                                   match_cost, mismatch_cost)
        idx = (score_a + score_b[::-1]).argmax()
        del score_a, score_b
        lhs = hirschberg(s[:mid], t[:idx], del_cost, ins_cost, match_cost,
                         mismatch_cost, vis, depth + 1)
        if vis:
            print("  " * depth, f"({s}, {t})")
        rhs = hirschberg(s[mid:], t[idx:], del_cost, ins_cost, match_cost,
                         mismatch_cost, vis, depth + 1)
        return [lhs[i] + rhs[i] for i in range(3)]


def test(s, t, del_cost, ins_cost, match_cost, mismatch_cost):
    score, *alignments = hirschberg(s, t, del_cost, ins_cost, match_cost,
                                    mismatch_cost, vis=False)
    true_score = align.globalmd(
        s, t,
        match_cost, mismatch_cost,
        ins_cost, ins_cost,
        del_cost, del_cost
    )[0].score
    assert score == true_score, (score, true_score)


def random_dna(length):
    return "".join(random.choices(("A", "G", "T", "C"), k=length))


def main():
    del_cost = -2
    ins_cost = -2
    match_cost = 2
    mismatch_cost = -1

    # rand tests
    n_tests = 100
    max_len = 20

    for i in range(n_tests):
        s = random_dna(random.randrange(max_len - 1) + 1)
        t = random_dna(random.randrange(max_len - 1) + 1)
        test(s, t, del_cost, ins_cost, match_cost, mismatch_cost)
    print("random tests passed")

    # vis test
    # left call is upper, right call is lower, offset represents depth

    n = 6
    a = random_dna(2 * n)
    b = random_dna(n)
    score, *alignments = hirschberg(a, b, del_cost, ins_cost, match_cost,
                                    mismatch_cost, vis=True)
    print(*["".join(s) for s in alignments], sep="\n")
    print(score)


if __name__ == "__main__":
    main()
