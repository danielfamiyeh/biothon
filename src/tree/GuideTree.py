from src.align.pairwise.AffineGlobal import *


class Node:
    """
    Node class, provides template for branching points in tree.
    """
    def __init__(self):
        """
        Node constructor.
        """
        self.weight = 1    # Taxon weight
        self.left = None   # Left child reference
        self.right = None  # Right child reference

    def set_weights(self, seqs):
        """
        Recursively assigns weights to nodes in all subtrees of node.
        Also assigns sequence to leaves connected to node.
        :param seqs:    List of sequences
        :return:        None
        """
        # If left child exists
        if self.left is not None:
            # Set its weight to current node's weight div 2
            self.left.weight = self.weight / 2

        # If the left child is another node
        if not isinstance(self.left, Leaf):
            # Perform recursive call to set_weights on left child
            self.left.set_weights(seqs)

        # Else, left child is a leaf in the tree
        else:
            # Assign a sequence to the leaf using
            #   its label as an index
            self.left.seq = seqs[int(self.left.label)]
            # Assign the left child a weight
            self.left.seq.weight = self.weight

        # If right child exists
        if self.right is not None:
            # Set its weight to the current node's weight div 2
            self.right.weight = self.weight / 2

        # If the right child is another node
        if not isinstance(self.right, Leaf):
            # Perform recursive call to set_weights on child
            self.right.set_weights(seqs)

        # Else, right child is a leaf in the tree
        else:
            # Assign the the child a sequence using
            #   the sequence's label as an index
            self.right.seq = seqs[int(self.right.label)]
            # Assign the right child a weight
            self.right.seq.weight = self.weight

    def postorder(self, func):
        """
        Performs post-order traversal of tree applying a function to every node.
        :param func:    Function to apply
        :return:        None
        """
        # Traverse right subtree
        try:
            self.right.postorder(func)
        # Perform function func on right leaf
        except AttributeError:
            func(self.right)

        # Traverse left subtree
        try:
            self.left.postorder(func)
        # Perform function func on left leaf
        except AttributeError:
            func(self.left)

    def preorder(self, func):
        """
        Performs pre-order traversal of tree applying a function to every node.
        :param func:    Function to apply
        :return:        None
        """
        # Traverse left subtree
        try:
            self.left.preorder(func)
        # Perform function func to left child
        except AttributeError:
            func(self.left)

        # Traverse right subtree
        try:
            self.right.preorder(func)
        # Perform function func on right child
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
        self.label = l  # Assign leaf a label => Denotes index in seq list
        self.weight = None
        self.seq = None # Leaf's sequence

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"({self.label} : {self.weight}) : {self.seq}"


class _Taxon:
    """
    Taxon class, represents sequences during guide tree construction.
    """
    def __init__(self, t_id, others):
        """
        Taxon constructor.
        :param t_id:    Numbered ID representing sequence's index in sequence list.
        :param others:  Number of other taxa in set of sequences.
        """
        self.id = t_id  # Taxon ID => Denotes index in seq list
        # Distance map used in cluster analysis
        self.dists = {str(i): -1 for i in range(others)}

        if t_id in self.dists:
            del self.dists[t_id] # Remove self from distance map

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
        :param mat:     Scoring matrix to define distance between sequences.
        """
        self.seqs = seqs    # Sequence set
        self.mat = mat      # Score matrix
        self.taxa = {str(i): _Taxon(str(i), len(seqs))
                     for i in range(len(seqs))}
        self.leaves = {str(i): Leaf(str(i))
                       for i in range(len(seqs))}
        self.nodes = {}

        self.root = None
        self.root_label = None
        self.start_label = None

        self.construct()

    def construct(self):
        """
        Constructs stringed, nested-tuple representation.
        :return: None
        """
        aligner = AffineGlobal(score_mat=self.mat)

        # Initialise tuple with highest distance = -1 and blank IDs
        #   highest[0]: distance = -1
        #   highest[1]: taxon #1 label = (blank)
        #   highest[2]: taxon #2 label = (blank)
        highest = -1, "", ""

        # Enumerate taxon map
        for i, t in enumerate(self.taxa):
            # Enumerate distances for current taxon
            for j, d in enumerate(self.taxa[t].dists):
                # Since distances are symmetric we only consider
                # entries above and including leading diagonal
                # (in a distance matrix)
                if int(d) > i:
                    # Calculate distance using NW
                    dist = aligner.align(self.seqs[int(t)], self.seqs[int(d)]).dist
                    # Assign distances to taxon pair
                    self.taxa[t].dists[d] = dist
                    self.taxa[d].dists[t] = dist
                    # Update highest tuple
                    if (i == 0 and j == 0) or dist < highest[0]:
                        highest = dist, t, d
        # Merge call to stringify first stage of cluster analysis
        # Also updates distance map and merges two taxa with
        #   the least distance
        self._merge(*highest)

        # Loops until cluster analysis is complete
        #   this is indicated by the taxa map having
        #   a length < 1 via successive merge calls
        while len(self.taxa) > 1:
            # Initialise highest distance tuple with
            # h[0]: distance = 0
            # h[1]: taxon #1 label = (blank)
            # h[2]: taxon #2 label = (blank)
            h = -1, "", ""

            # Enumerate over taxon left in taxon map
            for i, t in enumerate(self.taxa):
                # Enumerate over all distances of current taxon
                for j, other in enumerate(self.taxa[t].dists):
                    # Update highest
                    if (i == 0 and j == 0) or h[0] > self.taxa[t].dists[other]:
                        h = self.taxa[t].dists[other], t, other
            # Merge closest-related taxa
            #   and update taxon map
            self._merge(*h)

        # Tree rooting and weight assigning
        self.root = self.nodes[list(self.nodes.keys())[-1]]
        self.root.set_weights(self.seqs)
        root_str = str(self.root)

        # Get first index of label number
        first_digit = 0
        for i, c in enumerate(str(self.root)):
            if c.isdigit():
                first_digit = i
                break

        # Find colon separator to get index range of root label
        self.root_label = root_str[first_digit: root_str.find(":")-1]

        # Reverse tree string
        rev = "".join(reversed(str(self.nodes[list(self.nodes.keys())[-1]])))
        # Slice reversed tree string using colon separators
        rev = rev[rev.find(":") + 1:]
        rev = rev[rev.find(":") + 2:]
        # Get vicinity start label using sliced, reversed tree string
        self.start_label = "".join(reversed(rev[0:rev.find("(")]))
        # Slice start label using colon separator
        self.start_label = self.start_label[0:self.start_label.find(":")-1]

    def _merge(self, dist, i, j):
        """
        Merges two taxa in the distance map
        :param dist: Distance between two taxa
        :param i:   Destination taxon
        :param j:   Source taxon
        :return:    None
        """
        # New merged label using '/' notation from Oxford Academic
        new_id = f"({i}/{j})"
        l_depth = i.count("/")
        r_depth = j.count("/")
        depth = new_id.count("/")

        # Initialise new Node object
        node = Node()
        # Assign children nodes/leaves
        node.left = self.nodes[i] if l_depth > 0 else self.leaves[i]
        node.right = self.nodes[j] if r_depth > 0 else self.leaves[j]
        # Assign new ID to node
        self.nodes[new_id] = node

        # Delete pairs' distances
        del self.taxa[i].dists[j]
        del self.taxa[j].dists[i]

        # Average pairs' distances to other taxa
        values = [sum(x)/2 for x in zip(self.taxa[i].dists.values(),
                                        self.taxa[j].dists.values())]

        # Initialise new taxon based on merged pair
        new_taxon = _Taxon(new_id, len(values))
        # Assign distances to new taxon
        new_taxon.dists = dict(zip(self.taxa[i].dists.keys(), values))

        # Delete old taxa from taxon map
        del self.taxa[i]
        del self.taxa[j]

        # Insert new taxon into taxon map
        self.taxa[new_id] = new_taxon

        # Iterate over taxon map
        for t in self.taxa:
            # If current taxon is not newly merged taxon
            if t != new_id:
                # Then update distance to merged taxon
                self.taxa[t].dists[new_id] = (self.taxa[t].dists[i]
                                              + self.taxa[t].dists[j]) / 2
                # Delete old distances
                del self.taxa[t].dists[i]
                del self.taxa[t].dists[j]

