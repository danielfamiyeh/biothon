import unittest
from src.seq import *


class TestSeq(unittest.TestCase):
    def setUp(self):
        self.dna_sequence1 = "ATCCGCGTA"
        self.dna_sequence2 = "ATGGTCATA"
        self.rna_sequence = "AUGUCCAGA"
        self.prot_sequence = "MKLVMQV*"

        self.invalid_dna = "UAGGATCG"
        self.invalid_rna = "AGGCCAGT"
        self.invalid_prot = "MKLVMXQ*"

        self.dna1 = Seq(self.dna_sequence1)
        self.dna2 = Seq(self.dna_sequence2)
        self.rna = Seq(self.rna_sequence, SeqType.RNA)
        self.prot = Seq(self.prot_sequence, SeqType.PROTEIN)

    def testConstructor(self):
        self.assertEqual(str(self.dna1), self.dna_sequence1)
        self.assertEqual(str(self.dna2), self.dna_sequence2)
        self.assertEqual(str(self.rna), self.rna_sequence)
        self.assertEqual(str(self.prot), self.prot_sequence)

        self.assertEqual(self.dna1.seq_type, SeqType.DNA)
        self.assertEqual(self.rna.seq_type, SeqType.RNA)
        self.assertEqual(self.prot.seq_type, SeqType.PROTEIN)

    def testPointMutations(self):
        self.assertEqual(4, self.dna1.point_mutations(self.dna2))

    def testCount(self):
        count = self.dna1.count()
        self.assertEqual(2, count["A"])
        self.assertEqual(2, count["G"])
        self.assertEqual(3, count["C"])
        self.assertEqual(2, count["T"])

    def testComp(self):
        self.dna1.comp()
        self.rna.comp()

        self.assertEqual("TAGGCGCAT", str(self.dna1))
        self.assertEqual("UACAGGUCU", str(self.rna))
        self.assertEqual(self.dna_sequence1, str(~self.dna1))

        self.dna1.comp()
        self.rna.comp()

    def testRevComp(self):
        pass

    def testScribe(self):
        pass

    def testBackScribe(self):
        pass

    def testSlate(self):
        pass



if __name__ == '__main__':
    unittest.main()
