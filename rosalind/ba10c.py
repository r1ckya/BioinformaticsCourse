import sys

import numpy as np
from numpy.core.numeric import Infinity

input = sys.stdin.readline


def map_dict(d, f):
    return dict((k, f(v)) for k, v in d.items())


def viterbi(x, pi_states, transition_matrix, emission_matrix):

    n = len(x)
    m = len(pi_states)

    log_emission_matrix = map_dict(emission_matrix, np.log)
    log_transition_matrix = map_dict(transition_matrix, np.log)

    dp = np.zeros((n, m))
    logm = np.log(m)
    for i in range(m):
        dp[0, i] = -logm + log_emission_matrix[(pi_states[i], x[0])]

    for i in range(1, n):
        for j, cur_state in enumerate(pi_states):
            dp[i][j] = float("-inf")
            for k, prev_state in enumerate(pi_states):
                dp[i][j] = max(
                    dp[i][j],
                    dp[i - 1][k]
                    + log_transition_matrix[(prev_state, cur_state)]
                    + log_emission_matrix[(cur_state, x[i])],
                )

    states = [pi_states[np.argmax(dp[n - 1])]]
    for i in range(n - 2, -1, -1):
        cur_state = states[-1]
        k = np.argmax(
            [
                dp[i, k]
                + log_transition_matrix[(prev_state, cur_state)]
                + log_emission_matrix[(cur_state, x[i + 1])]
                for k, prev_state in enumerate(pi_states)
            ]
        )
        states.append(pi_states[k])
    return "".join(reversed(states))


def read():
    x = input().strip()
    input()

    x_states = input().split()
    input()

    pi_states = input().split()
    input()

    transition_matrix = dict()

    input()
    for s1 in pi_states:
        line = input().split()
        for s2, p in zip(pi_states, map(float, line[1:])):
            transition_matrix[(s1, s2)] = p
    input()

    emission_matrix = dict()

    input()
    for s1 in pi_states:
        line = input().split()
        for s2, p in zip(x_states, map(float, line[1:])):
            emission_matrix[(s1, s2)] = p

    return x, x_states, pi_states, transition_matrix, emission_matrix


def main():
    x, _, pi_states, transition_matrix, emission_matrix = read()
    print(viterbi(x, pi_states, transition_matrix, emission_matrix))


if __name__ == "__main__":
    main()
