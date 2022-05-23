import unittest

from src.matrix.ScoreMatrix import *
from src.parser import ScoreMatrixParser


class ScoreMatrixParserTest(unittest.TestCase):
    def test_something(self):
        pam250 = ScoreMatrixParser.load(AminoScoreType.PAM250)
        blosum62 = ScoreMatrixParser.load(AminoScoreType.BLOSUM62)

        for key in pam250:
            print(key, pam250[key], '\n')

        print('\n')

        for key in blosum62:
            print(key, blosum62[key], '\n')


if __name__ == '__main__':
    unittest.main()
