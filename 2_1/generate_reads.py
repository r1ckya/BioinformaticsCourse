import os
from enum import Enum
from functools import partial
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm


class Letters(Enum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"

    def to_tuple():
        return tuple(c.value for c in Letters)


def random_seq(vocab, length):
    return "".join(np.random.choice(vocab, size=length))


def random_dna(length):
    return random_seq(Letters.to_tuple(), length)


def get_counts(c, length):
    vocab = [x for x in Letters.to_tuple() if x != c]
    return dict(zip(vocab, np.random.multinomial(length, [1 / 3] * 3)))


def generate_read(start, length, dna, n_trials, e_low, e_high):
    res = []
    gt = []
    for i in range(length):
        j = (start + i) % len(dna)

        e = int(np.random.uniform(e_low, e_high))
        counts = get_counts(dna[j], e)
        counts[dna[j]] = n_trials - e
        c, cnt = max(counts.items(), key=lambda x: x[1])

        p = 1 - cnt / n_trials
        q = int(-10 * np.log10(p))

        res.append((c, chr(q + 33)))
        gt.append((dna[j], c == dna[j]))
    return res, gt


def generate_reads(
    dna,
    n_reads,
    n_trials,
    e_low,
    e_high,
    length_mean,
    length_std,
):
    lengths = np.random.normal(length_mean, length_std, n_reads).astype(int)
    starts = np.random.uniform(0, len(dna), size=n_reads).astype(int)

    with Pool(4) as p:
        res = p.starmap(
            partial(
                generate_read,
                dna=dna,
                n_trials=n_trials,
                e_low=e_low,
                e_high=e_high,
            ),
            tqdm(zip(starts, lengths), total=n_reads),
            chunksize=256,
        )
    return zip(*res)


def write_fastq(fname, reads):
    with open(fname, "w") as f:
        for i, read in enumerate(reads):
            print(f"@{i}", file=f)
            cs, qs = zip(*read)
            print(*cs, sep="", file=f)
            print("+", file=f)
            print(*qs, sep="", file=f)


def write_gt(fname, gts):
    with open(fname, "w") as f:
        for gt in gts:
            dna, correct = zip(*gt)
            print(*dna, sep="", file=f)
            print(*map(int, correct), sep="", file=f)


def main():
    dna_length = 50000
    n_reads = 50000
    n_trials = 100
    read_e_low = 1
    read_e_high = 76
    read_length_mean = 250
    read_length_std = 30

    dna = random_dna(dna_length)
    reads, gts = generate_reads(
        dna,
        n_reads,
        n_trials,
        e_low=read_e_low,
        e_high=read_e_high,
        length_mean=read_length_mean,
        length_std=read_length_std,
    )

    write_fastq("reads.fastq", reads)
    write_gt("gt.txt", gts)

    error_rate = []
    for gt in gts:
        _, correct = zip(*gt)
        error_rate.append(1 - sum(correct) / len(gt))

    mean_error_rate = sum(error_rate) / len(error_rate)

    print(f"mean_error_rate: {mean_error_rate}")


if __name__ == "__main__":
    main()
