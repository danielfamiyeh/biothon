import unittest
from src.matrix.ScoreMatrix import *
from src.align.pair.PairAligner import *

class TestPairAligner(unittest.TestCase):
    def setUp(self):
        self.seq1 = Seq("TGGTG")
        self.seq2 = Seq("ATCGT")
        self.aligner = PairAligner(score_mat=NucleoScoreMatrix(NucleoScoreType.IDENTITY))

    def testNeedlemanWunsch(self):
        self.assertEqual(3, self.aligner.needle(self.seq1, self.seq2).score)

    def testSmithWaterman(self):
        self.assertEqual(True, False)

    def testHirschbergs(self):
        self.assertEqual(True, False)

    def testBlast(self):
        self.assertEqual(True, False)

    def testFasta(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
