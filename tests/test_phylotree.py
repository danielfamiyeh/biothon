import unittest
from src.tree.PhyloTree import *
from src.markov.MarkovModel import *


class TestPhylogeneticTree(unittest.TestCase):
    def test_something(self):
        model = MarkovModel()
        self.sequences = [model.sequence(30) for _ in range(10)]
        for i, seq in enumerate(self.sequences):
            seq.id = str(i)
            print(repr(seq))

        tree = PhyloTree(self.sequences)
        tree.construct()
        print(tree)


if __name__ == '__main__':
    unittest.main()
