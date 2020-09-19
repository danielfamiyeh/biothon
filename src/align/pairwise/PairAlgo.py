from math import inf
from enum import Enum
from src.align.pairwise.RecurrenceRelation import *

class _Ptr(Enum):
    """
    Pointer for traceback algorithm
    """
    DIAG = 0
    UP = 1
    LEFT = 2
    NULL = 3


class PairAlgo:
    """
    Generalised pairwise alignment algorithm class.
    """
    def __init__(self, **kwargs):
        """
        :param kwargs:
            bound:  Set of 2 3-tuples (char, lambda, lambda)
            mat:    Set of matrix names
            rec:    Recurrence relations
            opt:    Method to get optimum value
            tb:     Traceback condition tuple
        """

        self.score_mat = kwargs.get("score_mat")

        # For affine gap penalty
        self.gap_open = kwargs.get("gap_open", -1)
        self.gap_extend = kwargs.get("gap_extend", 0)

        self.banded = kwargs.get("banded")
        self.boundary_conditions = kwargs.get("bound")
        self.matrix_names = kwargs.get("mat")
        self.r_relations = kwargs.get("rec")
        self.get_opt = kwargs.get("get_opt")
        self.tb_cond = kwargs.get("tb_cond")

    def align(self, s, t):
        rows = len(s) + 1
        cols = len(t) + 1
        scores = {"n": 0, "d": self.gap_open, "e": self.gap_extend}

        matrices = {}
        for i, name in enumerate(self.matrix_names):
            matrices[name] = [[-inf for _ in range(cols)]
                              for __ in range(rows)]

        tb_matrix = [[_Ptr.LEFT if i == 0 else
                      _Ptr.UP if j == 0 else
                      None for j in range(cols)]
                     for i in range(rows)]

        tb_matrix[0][0] = None
        matrices["T"] = tb_matrix

        for cond in self.boundary_conditions:
            mat = matrices[cond[0]]
            for i in range(rows):
                mat[i][0] = cond[1](i)
            for j in range(cols):
                mat[0][j] = cond[2](j)

        for i in range(1, rows):
            for j in range(1, cols):
                if abs(i-j) < (self.banded[0] + self.banded[1]):
                    scores["s"] = self.score_mat.lookup(s[i-1], t[j-1])

                    for key, rel_list in self.r_relations.items():
                        rr_values = []
                        for r in rel_list:
                            indices = r.get_indices(i, j)
                            if indices == -1:
                                continue
                            elif indices == 0:
                                rr_values.append(0)
                            else:
                                rr_values.append(matrices[r.maps_from][indices[0]][indices[1]]
                                                 + scores[r.score])

                        argmax = max(rr_values)
                        matrices[key][i][j] = argmax
                        if key == self.matrix_names[0]:
                            tb_matrix[i][j] = _Ptr(rr_values.index(argmax))

        for name in self.matrix_names:
            print(name)
            for row in matrices[name]:
                print(row)
            print("\n")

        for row in tb_matrix:
            print(row)

        i, j = self.get_opt(matrices[self.matrix_names[0]])
        aligned_s, aligned_t, pipes = "", "", ""

        char_count = 1

        while matrices[self.tb_cond[0]][i][j] != self.tb_cond[1]:
            current_ptr = tb_matrix[i][j]
            if current_ptr is _Ptr.DIAG:
                i -= 1
                j -= 1
                aligned_s += s[i]
                aligned_t += t[j]
                pipes += "|" if s[i] == t[j] else " "

            elif current_ptr is _Ptr.LEFT:
                j -= 1
                aligned_s += "-"
                aligned_t += t[j]
                pipes += " "

            else:
                i -= 1
                aligned_s += s[i]
                aligned_t += "-"
                pipes += " "

            if char_count % 60 == 0:
                aligned_s += "\n"
                aligned_t += "\n"
                pipes += "\n"
            char_count += 1

        print(aligned_s[::-1])
        print(pipes[::-1])
        print(aligned_t[::-1])