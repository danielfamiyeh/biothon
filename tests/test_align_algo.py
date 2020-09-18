import unittest

from src.align.pairwise.PairAlgo import *
from src.align.pairwise.Needle import *
from src.matrix.ScoreMatrix import *

class TestAlignAlgo(unittest.TestCase):
    def test_boundary_conditions(self):
        b_conditions = [("M", lambda i: - 2 * i,
                         lambda j: - 2 * j)]

        r_relations = (("M", -1, -1, "s"),
                       ("M", -1, 0, "d"),
                       ("M", 0, -1, "d"))

        get_opt = lambda matrix: (len(matrix)-1, len(matrix[0])-1)

        align_algo = PairAlgo(score_mat=NucleoScoreMatrix(NucleoScoreType.NON_UNIFORM),
                               mat=tuple("M"), bound=b_conditions, rec=r_relations, gap_open=-2,
                               get_opt=get_opt, tb_cond=("T", None))
        align_algo.align("AAAC", "AGC")

        needleAligner = Needle(score_mat=NucleoScoreMatrix(NucleoScoreType.NON_UNIFORM),
                               gap=-2)

        needleAligner.align("AAAC", "AGC")

if __name__ == '__main__':
    unittest.main()
