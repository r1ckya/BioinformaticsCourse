import sys

import numpy as np

input = sys.stdin.readline


def nanargmin2d(x):
    n = x.shape[0]
    idx = np.nanargmin(x)
    i = idx // n
    j = idx - i * n
    assert i < j
    return i, j


def main():
    n = int(input())
    dist = np.empty((2 * n - 2, 2 * n - 2))
    dist.fill(np.nan)

    for i in range(n):
        dist[i, :n] = list(map(int, input().split()))

    dist[np.diag_indices(n)] = np.nan

    ans = []
    denom = n - 2
    for m in range(n, 2 * n - 2):
        dist_q = (
            denom * dist
            - np.nansum(dist, axis=1, keepdims=True)
            - np.nansum(dist, axis=0, keepdims=True)
        )
        i, j = nanargmin2d(dist_q)
        delta = (np.nansum(dist[i]) - np.nansum(dist[j])) / denom
        ans.append((i, m, 0.5 * (dist[i, j] + delta)))
        ans.append((j, m, 0.5 * (dist[i, j] - delta)))
        dist[m, ...] = dist[..., m] = 0.5 * (
            dist[..., i] + dist[..., j] - dist[i, j]
        )
        denom -= 1
        dist[i, ...] = dist[..., i] = dist[j, ...] = dist[..., j] = np.nan

    i, j = nanargmin2d(dist)
    ans.append((i, j, dist[i, j]))

    ans.sort()
    for v, u, w in ans:
        print(f"{v}->{u}:{w:0.4f}")
    ans.sort(key=lambda x: x[1])
    for v, u, w in ans:
        print(f"{u}->{v}:{w:0.4f}")


if __name__ == "__main__":
    main()
