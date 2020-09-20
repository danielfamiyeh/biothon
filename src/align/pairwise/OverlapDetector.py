from src.align.pairwise.PairAlgo import *

def _get_opt(matrix):
    last_row = len(matrix)-1
    last_col = len(matrix[0])-1

    i_max, j_max, argmax = 0, 0, -inf
    for j in range(last_col + 1):
        if matrix[last_row][j] > argmax:
            i_max = last_row
            j_max = j
            argmax = matrix[i_max][j_max]

    for i in range(last_row + 1):
        if matrix[i][last_col] > argmax:
            i_max = i
            j_max = last_col
            argmax = matrix[i_max][j_max]

    return i_max, j_max


class OverlapDetector(PairAlgo):
    def __init__(self, **kwargs):
        gap = kwargs.get("gap", -1)

        # Recurrence Relations
        r_relations = {"M": [RecurrenceRelation(maps_from="M", i=-1, j=-1, score="s"),
                             RecurrenceRelation(maps_from="M", i=-1, j=0, score="d"),
                             RecurrenceRelation(maps_from="M", i=0, j=-1, score="d")]}

        # Boundary conditions
        b_conditions = [("M", lambda i: 0, lambda j: 0)]

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": tuple("M"),
                   "bound": b_conditions,
                   "banded": kwargs.get("banded", (inf, inf)),
                   "rec": r_relations,
                   "gap_open": gap,
                   "get_opt": _get_opt,
                   "tb_cond": ("T", None)}
        super().__init__(**_kwargs)
