import copy
import sys

from Bio import Phylo


def _fitch_up(clade):
    if clade.is_terminal():
        clade.name = {clade.name}
        return
    left_clade, right_clade = clade.clades
    _fitch_up(left_clade)
    _fitch_up(right_clade)

    intersection = left_clade.name & right_clade.name
    if len(intersection) > 0:
        clade.name = intersection
    else:
        clade.name = left_clade.name | right_clade.name


def _fitch_down(clade, c=None):
    if c is None or c not in clade.name:
        c = clade.name.pop()
    clade.name = c
    if clade.is_terminal():
        return
    left_clade, right_clade = clade.clades
    _fitch_down(left_clade, c)
    _fitch_down(right_clade, c)


def fitch(tree):
    tree = copy.deepcopy(tree)
    _fitch_up(tree.clade)
    _fitch_down(tree.clade)
    return tree


def main():
    tree = Phylo.read(sys.stdin, "newick")
    labeled_tree = fitch(tree)
    Phylo.draw(
        labeled_tree,
        label_colors=dict(A="#e6194b", C="#4363d8", G="#f58231", T="#911eb4"),
    )


if __name__ == "__main__":
    main()
