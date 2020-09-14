import unittest
from src.seq import *

class TestSeq(unittest.TestCase):
    def testConstructor(self):
        sequence = "ATTGCTGTCG"
        self.assertEqual(sequence, Seq(sequence))


if __name__ == '__main__':
    unittest.main()
