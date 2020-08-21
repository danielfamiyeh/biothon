from enum import Enum
from dist.Distance import *

alphabet_dna = "ATCG"      # DNA Nucleotides
alphabet_rna = "AUCG"      # RNA Nucleotides

# Single letter amino acid codes
alphabet_protein_1 = "ACDEFGHIKLMNPQRSTVWY*"

# Three letter amino acid codes
_alphabet_protein_3 = ["Ala", "Cys", "Asp", "Glu", "Phe", "Gly",
                       "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro",
                       "Gln", "Arg", "Ser", "Thr", "Val", "Trp",
                       "Tyr", "STO"]

# Full name amino acids
_alphabet_protein_full = ["Alanine", "Cysteine", "Aspartic Acid",
                          "Glutamic Acid", "Phenylalanine", "Glycine",
                          "Histidine", "Isoleucine", "Lysine", "Leucine",
                          "Methionine", "Asparagine", "Proline", "Glutamine",
                          "Arginine", "Serine", "Threonine", "Valine",
                          "Tryptophan", "Tyrosine", "STOP"]


# Amino acid lookup table
# Reads out index of amino acid in protein alphabet list
# E.g. To get index of start codon Methionine (AUG):
#   amino_acid_map["A"]["U"]["G"] reads out index 11
#   Methionine is stored at the 11th index in all protein alphabet lists
_amino_acid_map = {"U": {"U": {"U": 4, "C": 4, "A": 9, "G": 9},
                         "C": {"U": 15, "C": 15, "A": 15, "G": 15},
                         "A": {"U": 19, "C": 19, "A": -1, "G": -1},
                         "G": {"U": 1, "C": 1, "A": -1, "G": 18}},

                   "C": {"U": {"U": 9, "C": 9, "A": 9, "G": 9},
                         "C": {"U": 12, "C": 12, "A": 12, "G": 12},
                         "A": {"U": 6, "C": 6, "A": 13, "G": 13},
                         "G": {"U": 14, "C": 14, "A": 14, "G": 14}},

                   "A": {"U": {"U": 7, "C": 7, "A": 7, "G": 10},
                         "C": {"U": 16, "C": 16, "A": 16, "G": 16},
                         "A": {"U": 11, "C": 11, "A": 8, "G": 8},
                         "G": {"U": 15, "C": 15, "A": 14, "G": 14}},

                   "G": {"U": {"U": 17, "C": 17, "A": 17, "G": 17},
                         "C": {"U": 0, "C": 0, "A": 0, "G": 0},
                         "A": {"U": 2, "C": 2, "A": 3, "G": 3},
                         "G": {"U": 5, "C": 5, "A": 5, "G": 5}}}

class SeqType(Enum):
    """
    Biological sequence type enum
    """
    STRINGG = 0
    DNA = 1
    RNA = 2
    PROTEIN = 3


_nucleo_complement = {"A": "T", "C": "G", "T": "A", "G": "C"}
_nucleo_transition = {"A": "G", "C": "T", "G": "A", "T": "C"}

class Seq:
    """
    Sequence class, fundamental representation of biosequences with staistical
        and functional methods.
    """
    def __init__(self, seq, stype=SeqType.DNA, **kwargs):
        """
        Seq object constructor.
        :param seq:     Biosequence string.
        :param stype:   Sequence type
        :param kwargs:  name:   Name to identify Seq object (no spaces)
                        desc:   Description of sequence.
        """
        alpha = alphabet_dna if stype is SeqType.DNA else\
            alphabet_rna if stype is SeqType.RNA else\
            alphabet_protein_1

        self.seq = ""

        for i, c in enumerate(seq):
            if c.upper() not in alpha and c != "-":
                raise ValueError(f"Character {c} at index {i} of sequence is not valid for {stype}.")
            self.seq += c.upper()

        self.seq_type = stype
        self.name = kwargs.get("name", "").replace(" ", "")
        self.desc = kwargs.get("desc", "")
        self.label = ""                         # For Guide/PhyloTree creation - Do not use.
        self.weight = 1                         # For PhyloTree creation - Do not use.

    # Class Methods
    def transit(self, seq):
        """
        Transition-Transversion ratio method.
        Returns the transition/transversion ratio with respect to another DNA sequence of
            equal length.
        :param seq: Sequence to compare with.
        :return:    Ratio between transition subtitutions and transversion substituitons.
        """
        if len(self) == len(seq) and\
            self.seq_type is SeqType.DNA and\
                seq.seq_type is SeqType.DNA:
            transitions = 0
            transversions = 0
            for i in range(len(self)):
                if self.seq[i] != seq[i]:
                    if _nucleo_transition[self.seq[i]] == seq[i]:
                        transitions += 1
                    else:
                        transversions += 1
            return round(transitions/transversions, 11)
        else:
            raise ValueError("Sequences must both be DNA and be of equal"
                             "length for transition/transversion ratio.")

    def point_mutations(self, other_seq):
        """
        Counts the number of point mutations with respect to another sequence,
            of equal length and SeqType.
        :param other_seq:   Sequence to compare with
        :return:            Point mutation count.
        """
        if len(other_seq) == len(self):
            return Distance.hamming(self.seq, other_seq.seq)
        else:
            raise ValueError("Both sequences must be of same"
                             "length for point_mutations() method.")

    def find_motif(self, subseq):
        """
        Finds the indices where the subsequence given by subseq can be found.
        :param subseq:  Subsequence reprsenting motif.
        :return:        List of indices where motif can be found.
        """
        indices = []
        for i in range(len(self) - (len(subseq) - 1)):
            if self[i: i+len(subseq)] == subseq:
                indices.append(i)
        return indices

    def count(self):
        """
        Counts the number of residues in sequence.
        :return: Dictionary containing residue counts.
        """
        alpha = alphabet_dna if self.seq_type is SeqType.DNA else\
                alphabet_rna if self.seq_type is SeqType.RNA else SeqType.PROTEIN
        counter = {c: 0 for c in alpha}
        for c in self.seq:
            if c not in alpha:
                counter[c] = 0
            counter[c] += 1
        return counter

    def gc_content(self):
        """
        Calculates base composition (G-C count) of DNA sequence.
        :return:    G-C content.
        """
        if self.seq_type is SeqType.DNA:
            counter = self.count()
            gc = counter["G"] + counter["C"]
            return (gc / (gc + counter["A"] + counter["T"])) * 100
        else:
            raise TypeError("Sequence must be of DNA for GC content.")

    def comp(self):
        """
        DNA/RNA complement method
        :return: None
        """

        if self.seq_type is SeqType.DNA or self.seq_type is SeqType.RNA:
            self.seq = ''.join(["A" if base == "U" else _nucleo_complement[base]
                               for base in self.seq])
        else:
            raise TypeError("Sequence type must be DNA or RNA for complement "
                            "method.")

    def rev_comp(self):
        """
        DNA/RNA reverse complement method
        :return: None
        """

        if self.seq_type is SeqType.DNA or self.seq_type is SeqType.RNA:
            self.comp()
            self.seq = self.seq[::-1]
        else:
            raise TypeError("Sequence type must be DNA or RNA for reverse "
                            "complement method.")

    def scribe(self, *args):
        """
        DNA->RNA trranscirption method. Splices introns from sequence and converts,
            DNA sequence to (m)RNA sequence.
        :param args: Introns to be spliced.
        :return:     None
        """
        if self.seq_type is SeqType.DNA:
            self.seq_type = SeqType.RNA
            for s in args:
                self.seq = self.seq.replace(s, "")
            self.seq = ''.join(["U" if base == "T" else base for base in self.seq])
        else:
            raise TypeError("Sequence type must be DNA for transcription "
                            "method.")

    def back_scribe(self):
        """
        RNA->DNA back-transcription method.
        Converts RNA sequence back into DNA sequence.
        Spliced introns remain spliced.
        """
        if self.seq_type is SeqType.RNA:
            self.seq_type = SeqType.DNA
            self.seq = ''.join(["T" if base == "U" else base for base in self.seq])
        else:
            raise TypeError("Sequence type must be RNA for back-"
                            "transcription.")

    def slate(self, stop_codon=False, *args):
        """
        RNA translation method - turns RNA sequence into protein sequence.
        If called on a DNA sequence, sequence is transcribed first.

        :param stop_codon:  Method will halt after first stop codon if True
        :param *args:       Introns to be spliced if slate is called on a DNA
                                sequence.
        :return:  None
        """
        if self.seq_type is SeqType.DNA:
            self.scribe(*args)
            self.slate(stop_codon)

        elif self.seq_type is SeqType.RNA:
            codon = []
            chain = ""
            base_count = 0
            self.seq_type = SeqType.PROTEIN

            for base in self.seq:
                codon.append(base)
                base_count += 1

                if base_count % 3 == 0 and base_count > 0:
                    amino_index = _amino_acid_map[codon[0]][codon[1]][codon[2]]
                    if stop_codon and amino_index == -1:
                        break

                    chain += alphabet_protein_1[amino_index]
                    codon.clear()
            self.seq = chain

        else:
            raise TypeError(f"{self.seq_type} is not a translatable sequence "
                            f"type.")

    # Magic Methods
    def __repr__(self):
        string = f"{self.name} | {self.seq_type}:\n"
        string += "".join([f"{c}\n" if i > 0 and i % 60 == 0
                           else c for i, c in enumerate(self)])
        return string

    def __str__(self):
        return self.__repr__()

    def __invert__(self):
        if self.seq_type is SeqType.DNA or self.seq_type is SeqType.RNA:
            return Seq(''.join(["A" if base == "U" else _nucleo_complement[base]
                                for base in self.seq]), self.seq_type)
        else:
            raise TypeError("Sequence type must be DNA or RNA for complement "
                            "operator.")

    def __iter__(self):
        for char in self.seq:
            yield char

    def __add__(self, other):
        if self.seq_type is other.seq_type:
            return Seq(self.seq + other.seq, self.seq_type)
        else:
            raise TypeError("Sequences must be of the same type for "
                            "concatenation")

    def __iadd__(self, other):
        if type(self) is type(other):
            if self.seq_type is other.seq_type:
                self.seq += other.seq
                return self
            else:
                raise TypeError("Sequences must be of the same type for "
                                "(in-place) concatenation")
        else:
            raise TypeError("Param 'other' must be a Seq object too.")

    def __eq__(self, other):
        if type(self) is type(other):
            return self.seq == other.seq and self.seq_type is other.seq_type
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        try:
            return self.seq[key]
        except TypeError:
            return Seq(self.seq[key.start:key.stop], self.seq_type)
