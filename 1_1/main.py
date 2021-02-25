import glob
import random

from Bio import SeqIO
from distance import levenshtein, hamming


def hamming_distance(s, t):
    assert len(s) == len(t)
    res = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            res += 1
    return res


def hamming_search(s, t):  # search s in t
    n, m = len(s), len(t)
    assert n <= m
    ibest = 0
    dbest = m + 1
    for i in range(m - n + 1):
        cur = 0
        for j in range(n):
            if s[j] != t[i + j]:
                cur += 1
        if cur < dbest:
            dbest = cur
            ibest = i
    return ibest, t[ibest : ibest + n], dbest


def levenshtein_distance(s, t):
    n, m = len(s), len(t)
    if n > m:
        s, t = t, s
        n, m = m, n
    if n == 0:
        return m
    dp = [[i for i in range(n + 1)], [0] * (n + 1)]
    k = 0
    for i in range(m):
        k ^= 1
        dp[k][0] = i + 1
        for j in range(n):
            dp[k][j + 1] = min(
                dp[k ^ 1][j + 1] + 1,
                dp[k ^ 1][j] + (0 if t[i] == s[j] else 1),
                dp[k][j] + 1,
            )
    return dp[k][n]


def test(s, t, verbose=True):
    n = min(len(s), len(t))
    hdist = hamming_distance(s[:n], t[:n])
    assert hdist == hamming(s[:n], t[:n])
    levdist = levenshtein_distance(s, t)
    assert levdist == levenshtein(s, t)
    if verbose:
        print("hamming_distance:", hdist)
        print("hamming_search:", *hamming_search(s, t))
        print("levenshtein_distance:", levdist)
        print("-" * 100)


def random_dna(length):
    return "".join(random.choices(("A", "G", "T", "C"), k=length))


def main():
    # random tests
    n_tests = 100
    max_len = 100
    for i in range(n_tests):
        s = random_dna(random.randrange(max_len))
        t = random_dna(random.randrange(max_len))
        test(s, t, verbose=False)

    # seq tests
    for fname in glob.glob("data/*.fasta"):
        seqs = []
        for record in SeqIO.parse(fname, "fasta"):
            seqs.append(record.seq)
        print(fname)
        test(*seqs, verbose=True)


if __name__ == "__main__":
    main()
