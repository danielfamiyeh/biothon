import unittest

from src.seq.Kmer import *

class TestKmer(unittest.TestCase):
    def test_composition(self):
        comp3 = ["TAA", "AAT",
                "ATG", "TGC", "GCC",
                "CCA", "CAT", "ATG",
                "TGG", "GGG", "GGA",
                "GAT", "ATG", "TGT",
                "GTT"]
        string = "TAATGCCATGGGATGTT"
        self.assertEqual(comp3, get_composition(string, 3))


if __name__ == '__main__':
    unittest.main()
