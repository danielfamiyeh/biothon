import unittest

from src.seq import *
from src.align.pair.PairAlignment import *


class TestPairAlignment(unittest.TestCase):
    def setUp(self):
        self.aligned_seq1 = Seq("-TGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTG"
                                "TGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTGTGGTG", id=1)
        self.aligned_seq2 = Seq("ATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGT"
                                "ATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGT-", id=2345)
        self.score = 3
        self.alignment = PairAlignment(self.aligned_seq1,
                                       self.aligned_seq2, self.score)

    def testString(self):
        print(self.alignment)
        print(repr(self.alignment))

    def testScore(self):
        self.assertEqual(self.score, self.alignment.score)

    def testDistance(self):
        self.assertEqual(round(39.45, 3), self.alignment.dist)


if __name__ == '__main__':
    unittest.main()
