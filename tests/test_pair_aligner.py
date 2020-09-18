import unittest
from src.matrix.ScoreMatrix import *
from src.align.pair.PairAligner import *

class TestPairAligner(unittest.TestCase):
    def setUp(self):
        self.seq1 = Seq("TGGTG")
        self.seq2 = Seq("ATCGT")
        self.aligner = PairAligner(score_mat=NucleoScoreMatrix(NucleoScoreType.IDENTITY),
                                   match=1, mismatch=1, gap_open=3, gap_extend=1)

    def testNeedlemanWunsch(self):
        seq1 = "AAT"
        seq2 = "ACACT"
        string_seq = Seq(seq1, SeqType.STRING)
        string_seq2 = Seq(seq2, SeqType.STRING)
        print(self.aligner.needle_affine(string_seq, string_seq2))
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
