import unittest

from src.matrix.ScoreMatrix import *

class TestScoreMatrix(unittest.TestCase):
    def testBLOSUM62(self):
        matrix = AminoScoreMatrix(AminoScoreType.BLOSUM62)
        print(matrix)

    def testPAM250(self):
        matrix = AminoScoreMatrix(AminoScoreType.PAM250)
        print(matrix)

if __name__ == '__main__':
    unittest.main()
