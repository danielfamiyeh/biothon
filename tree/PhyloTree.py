from align.pair.PairAligner import PairAligner
from matrix.ScoreMatrix import *
from copy import deepcopy

class _Taxon:
    def __init__(self, t_id, others):
        """
        Taxon constructor
        :param t_id: Taxon ID
        :param others: Number of other taxa in tree.
        """
        self.id = t_id
        self.dists = {str(i): -1 for i in range(others)}
        if t_id in self.dists:
            del self.dists[t_id]

    def __repr__(self):
        return f"Taxon: {self.id}"

    def __str__(self):
        return self.__repr__()


class PhyloTree:
    def __init__(self, seqs):
        """
        Phylogenetic tree constructor
        :param seqs: Sequences to construct tree for
        """
        if len(seqs) > 3:
            self.S = deepcopy(seqs)
            self.taxa = {str(i) : _Taxon(str(i), len(seqs))
                         for i in range(len(seqs))}
            self.mergers = []
            self.tree = None
            self.trees = {}
        else:
            raise ValueError("Length of sequence list must be 3+"
                             "for tree construction.")

    def _construct(self):
        """
        Constructs nested tuple representation.
        :return: None
        """
        aligner = PairAligner(NucleoScoreMatrix(NucleoScoreType.IDENTITY))

        highest = -1, "", ""

        for i, t in enumerate(self.taxa):
            for j, d in enumerate(self.taxa[t].dists):
                if int(d) > i:
                    dist = aligner.needle(self.S[int(t)], self.S[int(d)]).dist
                    self.taxa[t].dists[d] = dist
                    self.taxa[d].dists[t] = dist
                    if (i == 0 and j == 0) or dist < highest[0]:
                        highest = dist, t, d
        self._merge(*highest)

        while len(self.taxa) > 1:
            h = -1, "", ""
            for i, t in enumerate(self.taxa):
                for j, other in enumerate(self.taxa[t].dists):
                    if (i == 0 and j == 0) or h[0] > self.taxa[t].dists[other]:
                        h = self.taxa[t].dists[other], t, other
            self._merge(*h)

        self.tree_repr = next(iter(self.taxa.keys()))

    def construct(self):
        """
        First constructs nested tuple representation of tree.
        Then this is parsed into actual tree structure.
        :return None
        """
        self._construct()
        self.tree = self.trees[self.mergers[len(self.mergers) - 1]]
        return self

    def _gen_leaf(self, n1, n2, d):
        return f"+{'-' * (2*d + (d-1)) + ' ' + str(n1)}\n" \
               f"│\n+{'-' * (2*d + (d-1)) + ' ' + str(n2)}"

    def _add_branch(self, *args, **kwargs):
        s = ""
        tree_str = self.trees[kwargs["subtree"]].split("\n")
        if "subtree2" not in kwargs:
            for i in range(len(tree_str)):
                if i == len(tree_str) - 2:
                    tree_str[i] = "+--" + tree_str[i]
                elif i == len(tree_str) - 1:
                    tree_str[i] = "│  " + tree_str[i]
                else:
                    tree_str[i] = " " * 3 + tree_str[i]

        if "singleton" in kwargs:
            sing_str = "\n+"
            sing_str += ("-" * (len(tree_str[0]) - 3)) + f" {kwargs['singleton']}"
            s += "\n".join(tree_str) + sing_str

        elif "n1" in kwargs and "n2" in kwargs:
            tuple2_str = ["\n+", "+"]
            tuple2_str[0] += ("-" * (len(tree_str[0]) - 3)) + f" {kwargs['n1']}"
            tuple2_str[1] += ("-" * (len(tree_str[0]) - 3)) + f" {kwargs['n2']}"
            s += "\n".join(tree_str) + "\n".join(tuple2_str)
        else:
            max1 = 0
            max2 = 0
            tree2_str = self.trees[kwargs["subtree2"]].split("\n")

            for line in tree_str:
                max1 = len(line) if len(line) > max1 else max1

            for line2 in tree2_str:
                max2 = len(line2) if len(line2) > max2 else max2

            diff = max1 - max2
            self._shift_subtree(tree2_str if diff > 0 else tree_str, abs(diff))

            last_row = tree_str[len(tree_str)-1].find("+")
            tree_str[len(tree_str)-1] = "-" * last_row + tree_str[len(tree_str)-1][last_row:]
            tree_str[len(tree_str)-1] = "+" + tree_str[len(tree_str)-1][1:len(tree_str[len(tree_str)-1])]

            min_branch = len(tree2_str[0])
            min_i = 0
            for i, line in enumerate(tree2_str):
                j = line.find("+")
                if j < min_branch:
                    min_i = i
                    min_branch = j

            for i in range(min_i):
                tree2_str[i] = "│" + tree2_str[i][1:]

            if min_branch != 0:
                tree2_str[min_i] = "+" + ("-" * (min_branch-1)) + tree2_str[min_i][min_branch:len(tree2_str[min_i])]

            s = "\n".join(tree_str) + "\n" + "\n".join(tree2_str)
        return s

    def _shift_subtree(self, st, n):
        """
        Returns a copy of the subtree shifted by n units.
        :param st: Subtree to shift
        :return: Shifted subtree
        """
        if n > 1:
            for i in range(len(st)):
                if i == len(st) - 2:
                    st[i] = f"{' ' * n}{st[i]}"
                elif i == len(st) - 1:
                    st[i] = f"{' ' * n}{st[i]}"
                else:
                    st[i] = f"{' ' * n}{st[i]}"


    def _parse(self, m):
        """
        Parses nested tuple representation of tree into printable tree structure
        using neighbour-joining history.
        :param m: Merger
        :return: None
        """

        for m in self.mergers:
            if self.trees[m] is None:
                # Singleton case
                if m[1].isdigit():
                    split_index = m.find("/")
                    self.trees[m] = self._add_branch(singleton=m[1:split_index],
                                                     subtree=m[split_index+1:len(m)-1])
                else:
                    # 2-tuple case
                    sp = m.split("/")
                    if m[2].isdigit() and sp[1][0].isdigit():
                        n1 = sp[0].split("(")[2]
                        n2 = sp[1][0:sp[1].find(")")]
                        subtree = m[m.find(")")+2:len(m)-1]

                        if self.trees[subtree] is not None:
                            self.trees[m] = self._add_branch(n1=n1, n2=n2,
                                                             subtree=subtree)

                    # Subtree case
                    else:
                        depth = 0
                        split_index = 0
                        for i in range(1, len(m)):
                            depth += 1 if m[i] == "(" else -1 if m[i] == ")" else 0
                            if depth == 0:
                                split_index = i + 1
                                break
                        self.trees[m] = self._add_branch(subtree=m[1:split_index],
                                                         subtree2=m[split_index+1:len(m)-1])

    def _merge(self, dist, i, j):
        """
        Merges two taxa in the distance map
        :param dist: Distance between two taxa
        :param i:   Destination taxon
        :param j:   Source taxon
        :return:    None
        """
        new_id = f"({i}/{j})"
        depth = new_id.count("/")
        self.trees[new_id] = self._gen_leaf(i, j, 1) if depth == 1 else None

        self.mergers.append(new_id)

        del self.taxa[i].dists[j]
        del self.taxa[j].dists[i]

        values = [sum(x)/2 for x in zip(self.taxa[i].dists.values(),
                                        self.taxa[j].dists.values())]

        new_taxon = _Taxon(new_id, len(values))
        new_taxon.dists = dict(zip(self.taxa[i].dists.keys(), values))

        del self.taxa[i]
        del self.taxa[j]

        self.taxa[new_id] = new_taxon

        for t in self.taxa:
            if t != new_id:
                self.taxa[t].dists[new_id] = (self.taxa[t].dists[i]
                                              + self.taxa[t].dists[j]) / 2
                del self.taxa[t].dists[i]
                del self.taxa[t].dists[j]
        self._parse(new_id)

    def __repr__(self):
        return ''.join([char + "\u0332" for char in "PhyloTree"]) +\
                f":\n{self.tree}"
