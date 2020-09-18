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
        r_relations = (("M", "M", -1, -1, "s"),
                       ("M", "M", -1, 0, "d"),
                       ("M", "M", 0, -1, "d"),
                       ("M", 0))

        # Boundary conditions
        b_conditions = [("M", lambda i: 0, lambda j: 0)]

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": tuple("M"),
                   "bound": b_conditions,
                   "rec": r_relations,
                   "gap_open": gap,
                   "get_opt": _get_opt,
                   "tb_cond": ("M", 0)}
        
        super().__init__(**_kwargs)
