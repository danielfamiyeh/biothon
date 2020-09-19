from src.align.pairwise.PairAlgo import *


def _get_opt(matrix):
    max_i, max_j, argmax = 0, 0, -inf
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > argmax:
                argmax = matrix[i][j]
                max_i = i
                max_j = j
    return max_i, max_j


class AffineLocal(PairAlgo):
    """
    Local aligner with affine gap penalty
    """
    def __init__(self, **kwargs):
        gap_open = kwargs.get("gap_open", -3)
        gap_extend = kwargs.get("gap_extend", -1)
        b_conditions = [("M", lambda i: 0,
                         lambda j: 0)]
        r_relations = {
            "F": [RecurrenceRelation(maps_from="M", i=-1, j=-1, score="s")],

            "G": [RecurrenceRelation(maps_from="M", i=-1, j=0, score="d"),
                  RecurrenceRelation(maps_from="G", i=-1, j=0, score="e")],

            "H": [RecurrenceRelation(maps_from="M", i=0, j=-1, score="d"),
                  RecurrenceRelation(maps_from="H", i=0, j=-1, score="e")],

            "M": [RecurrenceRelation(maps_from="F", i=0, j=0, score="n"),
                  RecurrenceRelation(maps_from="G", i=0, j=0, score="n"),
                  RecurrenceRelation(maps_from="H", i=0, j=0, score="n"),
                  RecurrenceRelation()]
        }

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": ("M", "F", "G", "H"),
                   "bound": b_conditions,
                   "banded": kwargs.get("banded", (inf, inf)),
                   "rec": r_relations,
                   "gap_open": gap_open,
                   "gap_extend": gap_extend,
                   "get_opt": _get_opt,
                   "tb_cond": ("M", 0)}

        super().__init__(**_kwargs)
