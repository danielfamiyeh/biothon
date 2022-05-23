import unittest

from src.dist.Distance import *

class TestDistance(unittest.TestCase):
    def setUp(self):
        self.str1 = "ATCGGTATCGA"
        self.str2 = "GATGGACGTGA"
        self.str3 = "ATACGA"
        self.hamming_distance = 7
        self.edit_distance1 = 6
        self.edit_distance2 = 5

    def testHamming(self):
        self.assertEqual(self.hamming_distance,
                         Distance.hamming(self.str1, self.str2))

    def testLevenshtein(self):
        self.assertEqual(self.edit_distance1,
                         Distance.edit(self.str1, self.str2))
        self.assertEqual(self.edit_distance2,
                         Distance.edit(self.str1, self.str3))


if __name__ == '__main__':
    unittest.main()
