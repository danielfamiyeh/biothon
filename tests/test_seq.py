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
        dna_rev_comp = (~self.dna1)[::-1]
        rna_rev_comp = (~self.rna)[::-1]

        self.dna1.rev_comp()
        self.rna.rev_comp()

        self.assertEqual(dna_rev_comp, str(self.dna1))
        self.assertEqual(rna_rev_comp, str(self.rna))

    def testScribe(self):
        self.dna1.scribe()
        self.dna2.scribe("GTC")
        transcribed = ''.join(["U" if base == "T" else base for base in self.dna1])
        spliced = "AUGAUA"

        print(self.dna1)

        self.assertEqual(SeqType.RNA, self.dna1.seq_type)
        self.assertEqual(SeqType.RNA, self.dna2.seq_type)

        self.assertEqual(transcribed, str(self.dna1))
        self.assertEqual(spliced, str(self.dna2))

    def testBackScribe(self):
        back_transcribed = ''.join(["T" if base == "U"
                                    else base for base in self.rna])
        self.rna.back_scribe()
        self.assertEqual(SeqType.DNA, self.rna.seq_type)
        self.assertEqual(back_transcribed, str(self.rna))

    def testSlate(self):
        translated = "MAMAPRTEINSTRING*"
        spliced = "MMPRTEISTRIG"
        longer_rna = Seq("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA",
                         SeqType.RNA)
        longer_dna = Seq("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA",
                         SeqType.RNA)
        longer_dna.back_scribe()

        longer_dna.slate(True, "GCC", "GCG", "AAT", "AAC")
        longer_rna.slate()

        self.assertEqual(SeqType.PROTEIN, longer_rna.seq_type)
        self.assertEqual(SeqType.PROTEIN, longer_dna.seq_type)

        self.assertEqual(translated, str(longer_rna))
        self.assertEqual(spliced, str(longer_dna))

    def testPrint(self):
        self.assertEqual(True, False)

    def testPrintFull(self):
        self.assertEqual(True, False)

if __name__ == '__main__':
    unittest.main()
