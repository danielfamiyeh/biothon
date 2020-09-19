from src.align.pairwise.PairAlgo import *
from src.align.pairwise.RecurrenceRelation import *


def _get_opt(matrix):
    max_i, max_j, argmax = 0, 0, -inf
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > argmax:
                argmax = matrix[i][j]
                max_i = i
                max_j = j
    return max_i, max_j

class LocalAligner(PairAlgo):
    def __init__(self, **kwargs):
        """
        Constructor
        :param kwargs:
            gap: (int) Gap penalty
            score_mat: (ScoreMatrix) Substitution matrix
        """

        gap = kwargs.get("gap")

        # Recurrence Relations
        r_relations = {
            "M": [RecurrenceRelation(maps_from="M", i=-1, j=-1, score="s"),
                  RecurrenceRelation(maps_from="M", i=-1, j=0, score="d"),
                  RecurrenceRelation(maps_from="M", i=0, j=-1, score="d"),
                  RecurrenceRelation()]
        }

        # Boundary conditions
        b_conditions = [("M", lambda i: 0, lambda j: 0)]

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": tuple("M"),
                   "bound": b_conditions,
                   "banded": kwargs.get("banded", (inf, inf)),
                   "rec": r_relations,
                   "gap_open": gap,
                   "get_opt": _get_opt,
                   "tb_cond": ("M", 0)}
        
        super().__init__(**_kwargs)
