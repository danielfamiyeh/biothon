from src.helper.RecRelation import *
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
        r_relations = (("M", "M", -1, -1, "s"),
                       ("M", "M", -1, 0, "d"),
                       ("M", "M", 0, -1, "d"),
                       ("M", 0))

        # Boundary conditions
        b_conditions = [("M", lambda i: gap * i,
                         lambda j: gap * j)]

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": tuple("M"),
                   "bound": b_conditions,
                   "rec": r_relations,
                   "gap_open": gap,
                   "get_opt": _get_opt,
                   "tb_cond": ("T", None)}
        super().__init__(**_kwargs)

    # Returns optimal score in alignment matrix

