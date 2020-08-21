from align.pair.PairAligner import *


class Node:
    """
    Node class, provides template for branching points in tree.
    as nodes in a tree.
    """
    def __init__(self):
        """
        Node constructor.
        """
        self.weight = 1
        self.left = None
        self.right = None

    def set_weights(self, seqs):
        """
        Recursively assigns weights to nodes in all subtrees of node.
        Also assigns sequence to leaves connected to node.
        :param seqs:    List of sequences
        :return:        None
        """
        if self.left is not None:
            self.left.weight = self.weight / 2

        if not isinstance(self.left, Leaf):
            self.left.set_weights(seqs)
        else:
            self.left.seq = seqs[int(self.left.label)]
            self.left.seq.weight = self.weight

        if self.right is not None:
            self.right.weight = self.weight / 2

        if not isinstance(self.right, Leaf):
            self.right.set_weights(seqs)
        else:
            self.right.seq = seqs[int(self.right.label)]
            self.right.seq.weight = self.weight

    def postorder(self, func):
        """
        Performs post-order traversal of tree applying a function to every node.
        :param func:    Function to apply
        :return:        None
        """
        try:
            self.right.postorder(func)
        except AttributeError:
            func(self.right)

        try:
            self.left.postorder(func)
        except AttributeError:
            func(self.left)

    def preorder(self, func):
        """
        Performs pre-order traversal of tree applying a function to every node.
        :param func:    Function to apply
        :return:        None
        """
        try:
            self.left.preorder(func)
        except AttributeError:
            func(self.left)

        try:
            self.right.preorder(func)
        except AttributeError:
            func(self.right)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"({self.left} : {self.right})"


class Leaf:
    """
    Leaf class, template for leaves in tree representing taxa/sequences.
    """
    def __init__(self, l):
        """
        Leaf object constructor.
        :param l: Leaf label number.
        """
        self.label= l
        self.weight = None
        self.seq = None

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"({self.label} : {self.weight}) : {self.seq}"


class _Taxon:
    """
    Taxon class, represents sequnces during guide tree construction.
    """
    def __init__(self, t_id, others):
        """
        Taxon constructor.
        :param t_id:    Numbered ID representing sequences index in sequence list.
        :param others:  Number of other taxa in set of sequences.
        """
        self.id = t_id
        self.dists = {str(i): -1 for i in range(others)}
        if t_id in self.dists:
            del self.dists[t_id]

    def __repr__(self):
        return f"Taxon: {self.id}"

    def __str__(self):
        return self.__repr__()


class GuideTree:
    """
    GuideTree class, used for phylogenetic approach in ClustalW algorithm.
    """
    def __init__(self, seqs, mat):
        """
        GuideTree constructor.
        :param seqs:    Set of sequences to construct tree for.
        :param mat:     Scoring matrix to define distance betweens sequences.
        """
        self.seqs = seqs
        self.mat = mat
        self.taxa = {str(i): _Taxon(str(i), len(seqs))
                     for i in range(len(seqs))}
        self.leaves = {str(i): Leaf(str(i)) for i in range(len(seqs))}
        self.nodes = {}
        self.root = None
        self.root_label = None
        self.start_label = None
        self.construct()

    def construct(self):
        """
        Constructs stringed. nested tuple representation.
        :return: None
        """
        aligner = PairAligner(self.mat)

        highest = -1, "", ""

        for i, t in enumerate(self.taxa):
            for j, d in enumerate(self.taxa[t].dists):
                if int(d) > i:
                    dist = aligner.needle(self.seqs[int(t)], self.seqs[int(d)]).dist
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

        self.root = self.nodes[list(self.nodes.keys())[-1]]
        self.root.set_weights(self.seqs)
        root_str = str(self.root)

        first_digit = 0
        for i, c in enumerate(str(self.root)):
            if c.isdigit():
                first_digit = i
                break

        self.root_label = root_str[first_digit: root_str.find(":")-1]

        rev = "".join(reversed(str(self.nodes[list(self.nodes.keys())[-1]])))
        rev = rev[rev.find(":") + 1:]
        rev = rev[rev.find(":") + 2:]

        self.start_label = "".join(reversed(rev[0:rev.find("(")]))
        self.start_label = self.start_label[0:self.start_label.find(":")-1]

    def _merge(self, dist, i, j):
        """
        Merges two taxa in the distance map
        :param dist: Distance between two taxa
        :param i:   Destination taxon
        :param j:   Source taxon
        :return:    None
        """
        new_id = f"({i}/{j})"
        l_depth = i.count("/")
        r_depth = j.count("/")
        depth = new_id.count("/")

        node = Node()
        node.left = self.nodes[i] if l_depth > 0 else self.leaves[i]
        node.right = self.nodes[j] if r_depth > 0 else self.leaves[j]
        self.nodes[new_id] = node

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

