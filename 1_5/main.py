from msa import multiple_sequence_alignment


def print_alignment(alignment):
    for a in alignment:
        print("".join(a))


def test(seqs, del_cost, ins_cost, match_cost, mismatch_cost):
    score, alignment = multiple_sequence_alignment(
        seqs, del_cost, ins_cost, match_cost, mismatch_cost
    )
    print("-" * 50)
    print(score)
    print_alignment(alignment)
    print("-" * 50)


def main():
    del_cost = -2
    ins_cost = -2
    match_cost = 2
    mismatch_cost = -1
    costs = (del_cost, ins_cost, match_cost, mismatch_cost)

    test(["ATGGA", "TGGAT", "GGATC"], *costs)
    test(["TTTAAA", "AAATTT", "GGGTTT"], *costs)
    test(["AAATTT", "GGGTTT", "AAAGGG"], *costs)
    test(["GAGAG", "AGAA", "GAGA"], *costs)
    test(["ACTGA", "GCATA", "GATGA"], *costs)
    test(["ACTA", "ACGTA"], *costs)
    test(["CTGA", "CATA", "ATGA", "ATAA"], *costs)


if __name__ == "__main__":
    main()
