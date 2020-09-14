from src.align.msa import *

class MultiAligner:
    def __init__(self, seqs, mat):
        """
        Multiple sequence aligner constructor.
        :param seqs: Set of sequences to align.
        :param mat: Scoring matrix for alignment.
        """
        if isinstance(seqs, list):
            if len(seqs) > 3:
                self.guide_tree = GuideTree(seqs, mat)
                self.tree_root = self.guide_tree.root
                self.profile = MSAProfile(self.guide_tree.leaves[
                                              self.guide_tree.start_label].seq, mat)
            else:
                raise ValueError("Param seqs must have three sequences or more.")
        else:
            raise TypeError("Param seqs must be of type list.")

    def _align(self, taxon):
        """
        Internal alignment method
        :param taxon: Taxon to align with profile
        :return: None
        """
        if taxon.label != self.guide_tree.start_label:
            self.profile.align(seq=taxon.seq)

    def align(self):
        """
        External alignment function that performs alignment to all
        taxa in guide tree in post-order.
        :return: None
        """
        self.tree_root.postorder(self._align)
