import sys
from collections import deque
from copy import deepcopy
from dataclasses import dataclass
from os.path import join

import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
from tqdm import tqdm


def window(s, n):
    for i in range(0, len(s) - n + 1):
        yield s[i : i + n]


def interpolate(a, b, c, d):
    return (a * b + c * d) / (b + d)


@dataclass
class EdgeInfo:
    coverage: float


class Edges:
    def __init__(self):
        self.edges = []
        self.deg_in = []
        self.deg_out = []

    def put(self, v, u, s, d=1):
        vd = self.edges[v]
        ud = vd.get(u)
        if ud is None:
            self.edges[v][u] = {s: EdgeInfo(coverage=d)}
            self.adjust_deg(v, u, +1)
        else:
            coverage = ud.get(s)
            if coverage is None:
                ud[s] = EdgeInfo(coverage=d)
                self.adjust_deg(v, u, +1)
            else:
                ud[s].coverage += d

    def remove(self, v, u, s):
        vd = self.edges[v]
        ud = vd.get(u)
        if ud is None:
            return
        ud.pop(s)
        self.adjust_deg(v, u, -1)
        if len(ud) == 0:
            vd.pop(u)

    def remove_all_s(self, v, u):
        vd = self.edges[v]
        if u in vd:
            self.adjust_deg(v, u, -len(vd.pop(u)))

    def adjust_deg(self, v, u, d):
        self.deg_out[v] += d
        self.deg_in[u] += d

    def add_vertex(self):
        self.edges.append(dict())
        self.deg_in.append(0)
        self.deg_out.append(0)

    def __getitem__(self, v):
        return self.edges[v]

    def __len__(self):
        return len(self.edges)

    def is_empty(self, v):
        return self.deg_in[v] == 0 and self.deg_out[v] == 0

    def remove_empty_vertexes(self):
        old_n_vertex = len(self.edges)
        old_edges = deepcopy(self.edges)

        delta = np.cumsum(
            [0] + [self.is_empty(v) for v in range(old_n_vertex)]
        )

        self.edges = []
        self.deg_in = []
        self.deg_out = []

        for _ in range(old_n_vertex - delta[-1]):
            self.add_vertex()

        for v in range(old_n_vertex):
            if delta[v + 1] - delta[v] == 0:
                for u, sd in old_edges[v].items():
                    for s, info in sd.items():
                        self.put(
                            v - delta[v + 1],
                            u - delta[u + 1],
                            s,
                            info.coverage,
                        )
        print(f"removed {delta[-1]} empty nodes", file=sys.stderr)


class DeBrujinGraph:
    def __init__(self, k):
        self.k = k
        self.edges = Edges()
        self.str2vertex = dict()

    def add_read(self, read):
        if len(read) < self.k:
            raise ValueError(f"found read shorter than {self.k}")

        for s in window(read, self.k):
            v = self._get_vertex(s[:-1])
            u = self._get_vertex(s[1:])
            self._add_edge(v, u, s)

    def _add_edge(self, v, u, s):
        self.edges.put(v, u, s)

    def _get_vertex(self, s):
        v = self.str2vertex.get(s)
        if v is None:
            v = len(self.str2vertex)
            self.str2vertex[s] = v
            self.edges.add_vertex()
        return v

    def compress_edge(self, v, u):
        if self.edges.deg_in[u] == 1 and self.edges.deg_out[u] == 1:
            ((to, sd),) = self.edges[u].items()
            ((uto_s, uto_info),) = sd.items()
            ((vu_s, vu_info),) = self.edges[v][u].items()
            self.edges.remove(v, u, vu_s)
            self.edges.remove(u, to, uto_s)
            self.edges.put(
                v,
                to,
                vu_s + uto_s[self.k - 1 :],
                interpolate(
                    vu_info.coverage,
                    len(vu_s),
                    uto_info.coverage,
                    len(uto_s),
                ),
            )
            return 1 + self.compress_edge(v, to)
        else:
            return 0

    def compress(self):
        compressed = 0
        for v in range(len(self.edges)):
            for u in tuple(self.edges[v].keys()):
                compressed += self.compress_edge(v, u)

        self.edges.remove_empty_vertexes()
        print(f"compressed total of {compressed} nodes", file=sys.stderr)


def remove_tails(dbg, remove_percentage=0.3):
    n_vertex = len(dbg.edges)
    scores = []

    for v in range(n_vertex):
        if dbg.edges.deg_in[v] > 0:
            for u, sd in dbg.edges[v].items():
                if dbg.edges.deg_in[u] == 1 and dbg.edges.deg_out[u] == 0:
                    ((s, info),) = sd.items()
                    scores.append((len(s) * info.coverage, v, u))

    scores.sort()
    print(len(scores))

    cut_idx = np.ceil(remove_percentage * len(scores)).astype(int)
    for _, v, u in scores[:cut_idx]:
        dbg.edges.remove_all_s(v, u)

    dbg.compress()


def is_connected(dbg, v, u):
    used = [False] * len(dbg.edges)
    queue = deque()
    used[v] = True

    while queue:
        v = queue.popleft()
        for to in dbg.edges[v]:
            if not used[to]:
                queue.append(to)
                used[to] = True
        if used[u]:
            return True

    return used[u]


def remove_bubbles(dbg):
    b = 0
    n_vertex = len(dbg.edges)
    for v in range(n_vertex):
        for u, sd in tuple(dbg.edges[v].items()):
            for s, info in tuple(sd.items()):
                if len(s) <= 2 * dbg.k:
                    dbg.edges.remove(v, u, s)
                    if not is_connected(dbg, v, u):
                        dbg.edges.put(v, u, s, info.coverage)
                    else:
                        b += 1
    dbg.compress()


def test(dbg):
    alls = []
    for v in range(len(dbg.edges)):
        for u, du in dbg.edges[v].items():
            for s, info in du.items():
                print(s, file=sys.stderr)
                alls.append(s)
                print(f"{v}->{u}", file=sys.stderr)
                print(f"coverage={info.coverage}", file=sys.stderr)
                print("-" * 100, file=sys.stderr)

    alls = "".join(map(str, alls))
    with open("dna.txt", "r") as f:
        dna = f.readline().strip()

    if alls == dna:
        print("assembly successful!")
    else:
        print("assembly failed!")


def plot_dbg(dbg, out_dir):
    scores = []
    n_vertex = len(dbg.edges)
    for v in range(n_vertex):
        for sd in dbg.edges[v].values():
            for info in sd.values():
                scores.append(info.coverage)

    with open(join(out_dir, "stats.txt"), "w") as f:
        print(f"number of nodes {len(dbg.edges)}", file=f)
        print(f"number of edges {len(scores)}", file=f)

    plt.hist(scores, range=[0, max(scores) + 1])
    plt.title("coverage distribution")
    plt.savefig(join(out_dir, "coverage.png"))
    plt.close()

    plt.hist(dbg.edges.deg_in)
    plt.title("deg in")
    plt.savefig(join(out_dir, "deg_in.png"))
    plt.close()

    plt.hist(dbg.edges.deg_out)
    plt.title("deg out")
    plt.savefig(join(out_dir, "deg_out.png"))
    plt.close()


def main():
    k = int(sys.argv[1])
    dbg = DeBrujinGraph(k)

    for record in tqdm(SeqIO.parse(sys.stdin, "fastq")):
        dbg.add_read(record.seq)

    plot_dbg(dbg, "original")

    dbg.compress()
    plot_dbg(dbg, "compressed")

    remove_tails(dbg)
    plot_dbg(dbg, "wo_tails")

    remove_bubbles(dbg)
    plot_dbg(dbg, "wo_bubbles")


if __name__ == "__main__":
    main()
