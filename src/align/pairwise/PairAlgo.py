from math import inf
from enum import Enum
from src.align.pairwise.PairAlign import *
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
            bound:  List of boundary conditions represented as 3-tuples
            mat:    Tuple of matrix names
            rec:    Dictionary of recurrence relation list for each matrix
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

        # Score added to recurrence relations
        scores = {"n": 0, "d": self.gap_open, "e": self.gap_extend}

        # Dictionary to get matrices by name
        matrices = {}

        # Initialise all matrices with negative infinity
        for i, name in enumerate(self.matrix_names):
            matrices[name] = [[-inf for _ in range(cols)]
                              for __ in range(rows)]

        # Traceback matrix
        tb_matrix = [[_Ptr.LEFT if i == 0 else
                      _Ptr.UP if j == 0 else
                      None for j in range(cols)]
                     for i in range(rows)]

        tb_matrix[0][0] = None
        matrices["T"] = tb_matrix

        # Initialise all matrices using boundary conditions
        for cond in self.boundary_conditions:
            mat = matrices[cond[0]]
            for i in range(rows):
                mat[i][0] = cond[1](i)
            for j in range(cols):
                mat[0][j] = cond[2](j)

        # DP algorithm
        for i in range(1, rows):
            for j in range(1, cols):
                # Conditional for banded alignment
                if abs(i-j) < (self.banded[0] + self.banded[1]):
                    scores["s"] = self.score_mat.lookup(s[i-1], t[j-1])

                    # Iterate over each tuple of recurrence relations
                    for key, rel_list in self.r_relations.items():
                        rr_values = []
                        for r in rel_list:
                            indices = r.get_indices(i, j)
                            # Conditional for local alignment case
                            if indices == 0:
                                rr_values.append(0)
                            else:
                                rr_values.append(matrices[r.maps_from][indices[0]][indices[1]]
                                                 + scores[r.score])

                        argmax = max(rr_values)
                        matrices[key][i][j] = argmax
                        if key == self.matrix_names[0]:
                            tb_matrix[i][j] = _Ptr(rr_values.index(argmax))

        # For debugging, uncomment to see alignment matrices and traceback matrix

 #       for name in self.matrix_names:
  #          print(name)
   #         for row in matrices[name]:
    #            print(row)
     #       print("\n")

      #  for row in tb_matrix:
       #     print(row)

        # Get indices of optimal score
        i, j = self.get_opt(matrices[self.matrix_names[0]])
        # Initialise alignment strings to blank strings
        aligned_s, aligned_t, pipes = [], [], []

        # Loop while traceback-halting condition has not been met
        while matrices[self.tb_cond[0]][i][j] != self.tb_cond[1]:
            current_ptr = tb_matrix[i][j]
            # Match/Mismatch case
            if current_ptr is _Ptr.DIAG:
                i -= 1
                j -= 1
                aligned_s += s[i]
                aligned_t += t[j]
                pipes += "|" if s[i] == t[j] else " "

            # Indel event in s-string
            elif current_ptr is _Ptr.LEFT:
                j -= 1
                aligned_s += "-"
                aligned_t += t[j]
                pipes += " "

            # Indel event in t-string
            else:
                i -= 1
                aligned_s += s[i]
                aligned_t += "-"
                pipes += " "

        aligned_s = aligned_s[::-1]
        pipes = pipes[::-1]
        aligned_t = aligned_t[::-1]

        for i in range(0, len(aligned_s), 60+i):
            if i > 0:
                aligned_s.insert(i, '\n')
                pipes.insert(i, '\n')
                aligned_t.insert(i, '\n')

        "".join(pipes)

        aligned_seq1 = s.copy()
        aligned_seq1.set_seq("".join(aligned_s))
        aligned_seq2 = t.copy()
        aligned_seq2.set_seq("".join(aligned_t))

        return PairAlign(aligned_seq1, aligned_seq2, pipes)