import unittest
from src.seq import *


class TestSeq(unittest.TestCase):
    def testConstructor(self):
        sequence = "ATTGCTGTCG"
        self.seq = Seq(sequence)
        self.assertEqual(sequence, self.seq.seq)


if __name__ == '__main__':
    unittest.main()
