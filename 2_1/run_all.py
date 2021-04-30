from generate_reads import main as generate_reads
from test_quake import main as test_quake
from trimmomatic import main as trimmomatic


def main():
    delim = "-" * 60
    print("generating reads...", flush=True)
    generate_reads()
    print(delim, flush=True)
    print("running trimmomatic...", flush=True)
    trimmomatic()
    print(delim, flush=True)
    print("running quake...", flush=True)
    test_quake()
    print(delim, flush=True)


if __name__ == "__main__":
    main()
