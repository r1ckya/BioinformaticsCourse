import sys

input = sys.stdin.readline


def read():
    seq = input().strip()

    input()
    states = input().split()
    input()
    input()

    transition_matrix = dict()

    for s1 in states:
        line = input().split()
        for s2, p in zip(states, map(float, line[1:])):
            transition_matrix[(s1, s2)] = p

    return seq, states, transition_matrix


def main():
    seq, states, transition_matrix = read()
    ans = 1 / len(states)
    for i in range(len(seq) - 1):
        ans *= transition_matrix[(seq[i], seq[i + 1])]
    print(ans)


if __name__ == "__main__":
    main()
