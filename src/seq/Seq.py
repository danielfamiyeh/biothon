from enum import Enum
from src.dist.Distance import *

alphabet_dna = "ATCG"  # DNA Nucleotides
alphabet_rna = "AUCG"  # RNA Nucleotides

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
    STRING = 0
    DNA = 1
    RNA = 2
    PROTEIN = 3


_dna_complement = {"A": "T", "C": "G", "T": "A", "G": "C"}
_rna_complement = {"A": "U", "C": "G", "U": "A", "G": "C"}
_dna_transition = {"A": "G", "C": "T", "G": "A", "T": "C"}


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

        # Initialise local alphabet based on SeqType
        alpha = alphabet_dna if stype is SeqType.DNA else \
            alphabet_rna if stype is SeqType.RNA else \
                alphabet_protein_1

        # Initialise sequence data to blank string
        self.seq = ""

        # Enumerate sequence string
        for i, c in enumerate(seq):
            # Check that character is part of alphabet
            if c.upper() not in alpha and c != "-":
                # Else throw an error
                raise ValueError(f"Character {c} at index {i} of sequence is not valid for {stype}.")
            # Append character to sequence string
            self.seq += c.upper()

        self.seq_type = stype
        self.id = str(kwargs.get("id", ""))
        self.name = kwargs.get("name", "")
        self.desc = kwargs.get("desc", "")
        self.label = ""  # For Guide/PhyloTree creation - Do not use.
        self.weight = 1  # For PhyloTree creation - Do not use.

    def transit(self, seq):
        """
        Transition-Transversion ratio method.
        Returns the transition/transversion ratio with respect to another DNA sequence of
            equal length.
        :param seq: Sequence to compare with.
        :return:    Ratio between transition subtitutions and transversion substituitons.
        """
        # Check that sequences are of the same length and type
        if len(self) == len(seq) and \
                self.seq_type is SeqType.DNA and \
                seq.seq_type is SeqType.DNA:
            # Initialise counters to zero
            transitions = 0
            transversions = 0
            # Iterate over sequence length
            for i in range(len(self)):
                # If characters do not match
                if self.seq[i] != seq[i]:
                    # Check if substitution is a transition
                    if _dna_transition[self.seq[i]] == seq[i]:
                        # If it is increment transition counter
                        transitions += 1
                    else:
                        # Else increment transversion counter
                        transversions += 1
            # Return rounded answer to 11 decimal places
            return round(transitions / transversions, 5)

        # If sequences are not compatible then throw error
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
        # Check that sequences are the same length
        if len(other_seq) == len(self):
            # Call hamming distance metric since this gives
            #   the number of point mutations
            return Distance.hamming(self.seq, other_seq.seq)
        # If sequences are not the same length then throw error
        else:
            raise ValueError("Both sequences must be of same"
                             "length for point_mutations() method.")

    def find_motif(self, subseq, **kwargs):
        """
        Finds the indices where the subsequence given by subseq can be found.
        :param subseq:  Subsequence reprsenting motif.
        :return:        List of indices where motif can be found.
        """
        overlap = kwargs.get("overlap", False)
        # Initialise empty list of indices
        indices = []
        # Iterate over sequence length wrt to substring length
        for i in range(0, len(self) - (len(subseq) - 1),
                       1 if overlap else len(subseq)):
            # If substring found
            if self[i: i + len(subseq)] == subseq:
                # Append index list
                indices.append(i)
        # Return index list
        return indices

    def count(self):
        """
        Counts the number of residues in sequence.
        :return: Dictionary containing residue counts.
        """
        # Initialise local alphabet variable based on SeqType
        alpha = alphabet_dna if self.seq_type is SeqType.DNA else \
            alphabet_rna if self.seq_type is SeqType.RNA else SeqType.PROTEIN

        # Initialise counter map
        counter = {c: 0 for c in alpha}

        # Iterate over characters in sequence
        for c in self.seq:
            if c not in alpha:
                counter[c] = 0
            # Increment counter map
            counter[c] += 1
        # Return counter map
        return counter

    def gc_content(self):
        """
        Calculates base composition (G-C count) of DNA sequence.
        :return:    G-C content.
        """
        if self.seq_type is SeqType.DNA:
            # Count all residues
            counter = self.count()
            # Sum GC counts
            gc = counter["G"] + counter["C"]
            # Return GC content
            return (gc / (gc + counter["A"] + counter["T"])) * 100
        else:
            # Throw type error if sequence isn't a DNA sequence
            raise TypeError("Sequence must be of DNA for GC content.")

    def comp(self):
        """
        DNA/RNA complement method
        :return: None
        """

        # Check that sequence is DNA or RNA
        if self.seq_type is SeqType.DNA:
            self.seq = ''.join([_dna_complement[base] for base in self.seq])
        elif self.seq_type is SeqType.RNA:
            self.seq = ''.join(_rna_complement[base] for base in self.seq)
        else:
            # Throw type error if sequence is not of valid type
            raise TypeError("Sequence type must be DNA or RNA for complement "
                            "method.")

    def rev_comp(self):
        """
        DNA/RNA reverse complement method
        :return: None
        """

        # Check that sequence is DNA or RNA
        if self.seq_type is SeqType.DNA or self.seq_type is SeqType.RNA:
            # Get complement
            self.comp()
            # Reverse complement
            self.seq = self.seq[::-1]
        else:
            # Throw type error if sequence is not of valid type
            raise TypeError("Sequence type must be DNA or RNA for reverse "
                            "complement method.")

    def scribe(self, *args):
        """
        DNA->RNA trranscirption method. Splices introns from sequence and converts,
            DNA sequence to (m)RNA sequence.
        :param args: Introns to be spliced.
        :return:     None
        """
        # Check that sequence is a DNA sequence
        if self.seq_type is SeqType.DNA:
            # Change sequence to RNA
            self.seq_type = SeqType.RNA
            if len(args) > 0:
                spliced = ""
                for i in range(0, len(self)-2, 3):
                    codon = self[i:i+3]
                    spliced += "" if codon in args else codon

                self.seq = spliced
            # Join spliced sequence
            self.seq = ''.join(["U" if base == "T" else base for base in self.seq])
        else:
            # Throw error if sequence is not a DNA sequence
            raise TypeError("Sequence type must be DNA for transcription "
                            "method.")

    def back_scribe(self):
        """
        RNA->DNA back-transcription method.
        Converts RNA sequence back into DNA sequence.
        Spliced introns remain spliced.
        """
        # Check that sequence is an RNA sequence
        if self.seq_type is SeqType.RNA:
            # Change sequence type back to DNA
            self.seq_type = SeqType.DNA
            # Perform back transcription, introns are lost
            self.seq = ''.join(["T" if base == "U" else base for base in self.seq])
        else:
            # Throw error if sequence is not an RNA sequence
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
        # If sequence is DNA sequence
        if self.seq_type is SeqType.DNA:
            # Transcribe first, then translate via recursive call
            self.scribe(*args)
            self.slate(stop_codon)

        # If sequence is RNA
        elif self.seq_type is SeqType.RNA:
            # Initialise empty codon list
            codon = []
            # Initialise blank chain string
            chain = ""
            # Initialse base (residue) count to zero
            base_count = 0
            # Change sequence type to protein
            self.seq_type = SeqType.PROTEIN

            # Iterate over residues in sequence
            for base in self.seq:
                # Append current base to codon list
                codon.append(base)
                base_count += 1

                # If we have three residues in codon list
                if base_count % 3 == 0 and base_count > 0:
                    # Lookup amino acid using residues and store
                    #   in amino_index variable
                    amino_index = _amino_acid_map[codon[0]][codon[1]][codon[2]]
                    if stop_codon and amino_index == -1:
                        break

                    # Add amino acid to polypeptide chain
                    chain += alphabet_protein_1[amino_index]
                    # Clear codon list
                    codon.clear()
            self.seq = chain

        else:
            # Throw error if sequence is not translatable
            raise TypeError(f"{self.seq_type} is not a translatable sequence "
                            f"type.")

    def copy(self):
        return Seq(self.seq, self.seq_type, id=self.id,
                   name=self.name, desc=self.desc)

    def set_seq(self, seq):
        self.seq = seq

    def __repr__(self):
        string = f"Seq({str(self)})"
        string += "\n"
        string += f"{self.id + '|' if len(self.id) > 0 else ''}" \
                 f"{self.name + '|' if len(self.name) > 0 else ''}" \
                 f"{self.desc if len(self.desc) > 0 else ''}"

        return string

    def __str__(self):
        return "".join([f"{c}\n" if i > 0 and i % 60 == 0
                           else c for i, c in enumerate(self)])
    def __invert__(self):
        if self.seq_type is SeqType.DNA:
            return Seq(''.join([_dna_complement[base] for base in self]),
                       self.seq_type, id=self.id, name=self.name, desc=self.desc)

        elif self.seq_type is SeqType.RNA:
            return Seq(''.join([_rna_complement[base] for base in self]),
                       self.seq_type, id=self.id, name=self.name, desc=self.desc)
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
