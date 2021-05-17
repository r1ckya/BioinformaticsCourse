import sys
from dataclasses import dataclass

from Bio import SeqIO


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

        print(f"compressed total of {compressed} nodes", file=sys.stderr)


def main():
    k = int(sys.argv[1])
    dbg = DeBrujinGraph(k)

    for record in SeqIO.parse("reads.fastq", "fastq"):
        dbg.add_read(record.seq)

    dbg.compress()

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


if __name__ == "__main__":
    main()
