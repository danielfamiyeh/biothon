import unittest
from src.markov.MarkovModel import *


class TestMarkov(unittest.TestCase):
    def setUp(self):
        self.tm = [[0.423, 0.151, 0.168, 0.258],
                   [0.399, 0.184, 0.063, 0.354],
                   [0.314, 0.189, 0.176, 0.321],
                   [0.258, 0.138, 0.187, 0.415]]
        self.pi = [0.345, 0.158, 0.159, 0.337]

    def testWithNoParams(self):
        model = MarkovModel()
        seq = model.sequence()
        self.assertEqual(Seq, type(seq))
        self.assertEqual(SeqType.DNA, seq.seq_type)
        self.assertTrue(1 <= len(seq) <= 30)

    def testModelType(self):
        dna_model = MarkovModel()
        rna_model = MarkovModel(type=SeqType.RNA)
        protein_model = MarkovModel(type=SeqType.PROTEIN)

        dna_seq = dna_model.sequence()
        rna_seq = rna_model.sequence()
        protein_seq = protein_model.sequence()

        self.assertEqual(SeqType.DNA, dna_seq.seq_type)
        self.assertEqual(SeqType.RNA, rna_seq.seq_type)
        self.assertEqual(SeqType.PROTEIN, protein_seq.seq_type)

    def testWithOnlyTransitionMatrix(self):
        model = MarkovModel(tm=self.tm)
        for i, row in enumerate(model.tm):
            self.assertEqual(self.tm[i], row)
        model.sequence()

    def testWithSteadyStateVector(self):
        model = MarkovModel(pi=self.pi)
        self.assertEqual(self.pi, model.pi)
        model.sequence()

    def testSequenceLength(self):
        model = MarkovModel(type=SeqType.RNA,
                            pi=self.pi,
                            tm=self.tm)
        seq = model.sequence(40)
        self.assertEqual(40, len(seq))


if __name__ == '__main__':
    unittest.main()
