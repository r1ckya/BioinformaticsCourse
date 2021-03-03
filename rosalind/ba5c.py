import sys

input = sys.stdin.readline


def main():
    s = input().strip()
    t = input().strip()

    n, m = len(s), len(t)
    dp = [[0] * (m + 1) for i in range(n + 1)]
    for i in range(n):
        for j in range(m):
            if s[i] == t[j]:
                dp[i + 1][j + 1] = dp[i][j] + 1
            else:
                dp[i + 1][j + 1] = max(dp[i + 1][j], dp[i][j + 1])

    cur_i = n - 1
    cur_j = m - 1
    ans = []
    while cur_i >= 0 and cur_j >= 0:
        if s[cur_i] == t[cur_j]:
            ans.append(s[cur_i])
            cur_i -= 1
            cur_j -= 1
        elif dp[cur_i + 1][cur_j + 1] == dp[cur_i][cur_j + 1]:
            cur_i -= 1
        else:
            cur_j -= 1

    print("".join(reversed(ans)))


if __name__ == "__main__":
    main()
