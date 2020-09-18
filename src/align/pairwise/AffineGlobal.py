from src.align.pairwise.PairAlgo import *


def _get_opt(matrix):
    return len(matrix)-1, len(matrix[0])-1


class AffineGlobal(PairAlgo):
    """
    Global aligner with affine gap penalty
    """
    def __init__(self, **kwargs):
        gap_open = kwargs.get("gap_open", -3)
        gap_extend = kwargs.get("gap_extend", -1)
        b_conditions = [("M", lambda i: gap_open + i*gap_extend,
                         lambda j: gap_open + j*gap_extend)]

        r_relations = {
            "F": [("M", -1, -1, "s")],
            "G": [("M", -1, 0, "d"), ("G", -1, 0, "e")],
            "H": [("M", 0, -1, "d"), ("H", 0, -1, "e")],
            "M": [("F", 0, 0, "n"), ("G", 0, 0, "n"), ("H", 0, 0, "n")]
        }

        _kwargs = {"score_mat": kwargs.get("score_mat"),
                   "mat": ("M", "F", "G", "H"),
                   "bound": b_conditions,
                   "rec": r_relations,
                   "gap_open": gap_open,
                   "gap_extend": gap_extend,
                   "get_opt": _get_opt,
                   "tb_cond": ("T", None)}

        super().__init__(**_kwargs)
