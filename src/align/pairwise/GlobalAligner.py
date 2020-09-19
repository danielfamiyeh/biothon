from src.align.pairwise.PairAlgo import *


def _get_opt(matrix):
    return len(matrix)-1, len(matrix[0])-1


class GlobalAligner(PairAlgo):
    """
    Needleman-Wunsch aligner
    """
    def __init__(self, **kwargs):
        """
        :param kwargs:
            gap: (int) Gap penalty
            score_mat: (ScoreMatrix) Substitution matrix
        """
        gap = kwargs.get("gap")

        # Recurrence Relations
        r_relations = {"M": [RecurrenceRelation(maps_from="M", i=-1, j=-1, score="s"),
                             RecurrenceRelation(maps_from="M", i=-1, j=0, score="d"),
                             RecurrenceRelation(maps_from="M", i=0, j=-1, score="d")]}

        # Boundary conditions
        b_conditions = [("M", lambda i: gap * i,
                         lambda j: gap * j)]

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": tuple("M"),
                   "bound": b_conditions,
                   "banded": kwargs.get("banded", (inf, inf)),
                   "rec": r_relations,
                   "gap_open": gap,
                   "get_opt": _get_opt,
                   "tb_cond": ("T", None)}
        super().__init__(**_kwargs)

    # Returns optimal score in alignment matrix

