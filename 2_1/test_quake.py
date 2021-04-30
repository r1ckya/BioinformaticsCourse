import os
from collections import defaultdict


def read_gts(fname):
    with open(fname, "r") as f:
        nts, correct = [], []
        for i, line in enumerate(f):
            if i % 2 == 0:
                nts.append(line.strip())
            else:
                correct.append(tuple(map(int, line.strip())))
        return nts, correct


def read_correction_log(fname):
    with open(fname, "r") as f:
        res = defaultdict(list)
        idx = 0
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                idx = int(line[1:])
            else:
                _, pos, new_nt, old_nt = line.split()
                res[idx].append((int(pos) - 1, new_nt, old_nt))
    return res


def run_quake():
    kmer_size = 12
    reads_file = "reads.fastq"
    counts_file = "counts.txt"
    cutoff_threshold = 0.0002

    os.system(
        f"""Quake/bin/count-qmers -k {kmer_size}"""
        f""" -f {reads_file} > {counts_file}"""
    )
    os.system(
        f"""Quake/bin/correct -r {reads_file} -k {kmer_size}"""
        f""" -c {cutoff_threshold} -m {counts_file} --log -p 4"""
    )


def main():
    run_quake()

    nts, correct = read_gts("gt.txt")
    correction_log = read_correction_log("reads.fastq.log")

    correct_corrections = 0
    wrong_nt_corrections = 0
    wrong_pos_corrections = 0

    for idx, corrections in correction_log.items():
        for pos, new_nt, old_nt in corrections:
            if correct[idx][pos] == 0:
                if nts[idx][pos] == new_nt:
                    correct_corrections += 1
                else:
                    wrong_nt_corrections += 1
            else:
                wrong_pos_corrections += 1

    print(f"correct_corrections: {correct_corrections}")
    print(f"wrong_nt_corrections: {wrong_nt_corrections}")
    print(f"wrong_pos_corrections: {wrong_pos_corrections}")


if __name__ == "__main__":
    main()
