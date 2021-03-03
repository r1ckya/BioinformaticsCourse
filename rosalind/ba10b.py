import sys

input = sys.stdin.readline


def read():
    x = input().strip()

    input()
    x_states = input().split()
    input()

    pi = input().strip()
    input()
    pi_states = input().split()
    input()

    input()

    emission_matrix = dict()

    for s1 in pi_states:
        line = input().split()
        for s2, x in zip(x_states, map(float, line[1:])):
            emission_matrix[(s1, s2)] = x

    return x, pi, emission_matrix


def main():
    x, pi, emission_matrix = read()
    ans = 1
    for e_x, e_pi in zip(x, pi):
        ans *= emission_matrix[(e_pi, e_x)]
    print(ans)


if __name__ == "__main__":
    main()
