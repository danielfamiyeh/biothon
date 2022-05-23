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

    def testGCContent(self):
        self.assertEqual(100 * (5/9), self.dna1.gc_content())

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

    def testFindMotif(self):
        seq = Seq("AAAAAA")
        self.assertEqual([0, 1, 2, 3], seq.find_motif("AAA", overlap=True))
        self.assertEqual([0, 3], seq.find_motif("AAA"))

    def testTransitTransver(self):
        self.assertEqual(round(1/3, 5), self.dna1.transit(self.dna2))

    def testPrint(self):
        self.dna1.id = "1234634"
        self.dna1.name = "Name"
        self.dna1.desc = "Short description of sequence."
        repr_split = repr(self.dna1).split("\n")

        self.assertEqual(str(self.dna1), self.dna_sequence1)
        self.assertEqual("Seq("+self.dna_sequence1+")", repr_split[0])
        self.assertEqual("1234634|Name|Short description of sequence.", repr_split[1])

    def testSplit(self):
        split_seqs = self.dna1.split("G")
        self.assertEqual("ATCC", split_seqs[0].seq)
        self.assertEqual("C", split_seqs[1].seq)
        self.assertEqual("TA", split_seqs[2].seq)

    def testMagicMethods(self):
        # Iter test
        for i, char in enumerate(self.dna1):
            self.assertEqual(char, self.dna_sequence1[i])

        # Concat test
        self.assertEqual(self.dna_sequence1 + self.dna_sequence2,
                         str(self.dna1 + self.dna2))

        # Inplace concat test
        self.dna1 += self.dna2
        self.assertEqual(self.dna_sequence1 + self.dna_sequence2,
                         str(self.dna1))

        # Equality tests
        self.assertEqual(self.dna1, self.dna1)
        self.assertNotEqual(self.dna1, self.dna2)

        # Length test
        self.assertEqual(len(self.dna_sequence2), len(self.dna2))

        # Getitem tests
        self.assertEqual(self.dna_sequence1[4], self.dna1[4])
        self.assertEqual(self.dna_sequence2[1:3], self.dna2[1:3])

if __name__ == '__main__':
    unittest.main()
