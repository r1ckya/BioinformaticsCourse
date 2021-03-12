import sys

import numpy as np

input = sys.stdin.readline


def normalize_rows(x):
    return x / x.sum(axis=1, keepdims=True)


def read():
    x = input().strip()
    input()

    x_states = input().split()

    input()
    y = input().strip()
    input()

    y_states = input().split()

    return x, x_states, y, y_states


def round3(x):
    return (x * 1e3 + 0.5).astype(int) / 1e3


def main():
    x, x_states, y, y_states = read()

    transition_matrix = np.zeros((len(y_states), len(y_states)))
    emission_matrix = np.zeros((len(y_states), len(x_states)))

    x_state_to_idx = dict((k, v) for v, k in enumerate(x_states))
    y_state_to_idx = dict((k, v) for v, k in enumerate(y_states))

    for a, b in zip(x, y):
        i = x_state_to_idx[a]
        j = y_state_to_idx[b]
        emission_matrix[j, i] += 1

    for i in range(len(y) - 1):
        a = y_state_to_idx[y[i]]
        b = y_state_to_idx[y[i + 1]]
        transition_matrix[a, b] += 1

    transition_matrix[transition_matrix.sum(1) == 0, 0] = 1
    emission_matrix[emission_matrix.sum(1) == 0, 0] = 1
    transition_matrix = normalize_rows(transition_matrix)
    emission_matrix = normalize_rows(emission_matrix)

    print(" \t" + "\t".join(y_states))
    for s, row in zip(y_states, transition_matrix):
        print(f"{s}\t", end="")
        print(*round3(row), sep="\t")
    print("--------")
    print(" \t" + "\t".join(x_states))
    for s, row in zip(y_states, emission_matrix):
        print(f"{s}\t", end="")
        print(*round3(row), sep="\t")


if __name__ == "__main__":
    main()
