import os


def read_gts(fname):
    with open(fname, "r") as f:
        return [
            tuple(map(int, line.strip()))
            for i, line in enumerate(f)
            if i % 2 == 1
        ]


def main():
    os.system(
        f"""java -jar trimmomatic-0.39.jar SE -phred33 -quiet"""
        """ -trimlog trimlog.txt reads.fastq outreads.fastq"""
        """ LEADING:3 TRAILING:3 SLIDINGWINDOW:10:3 MINLEN:36"""
    )

    gts = read_gts("gt.txt")

    with open("trimlog.txt", "r") as f:
        trimlog = [tuple(map(int, line.split())) for line in f]

    tp = 0
    fp = 0
    removed = 0
    trimmed = 0
    validated = 0

    for log, gt in zip(trimlog, gts):
        _, length, start, end, tail_length = log
        if length == 0:
            n_true = sum(gt)
            fp += n_true
            tp += len(gt) - n_true
            removed += 1
        else:
            if length == len(gt):
                validated += 1
            else:
                trimmed += 1

            n_true = sum(gt[:start])
            fp += n_true
            tp += start - n_true

            n_true = sum(gt[end:])
            fp += n_true
            tp += tail_length - n_true

    print(f"tp_ratio: {tp / (tp + fp)}")  # ~ mean error ratio
    print(f"fp_ratio: {fp / (tp + fp)}")

    print(f"validated: {validated}")
    print(f"trimmed: {trimmed}")
    print(f"removed: {removed}")

    # choose thr 0.0005


if __name__ == "__main__":
    main()
