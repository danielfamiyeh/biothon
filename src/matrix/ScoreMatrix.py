from enum import Enum
from src.parser.ScoreMatrixParser import *


class _ScoreMatrix:
    """
    Inner ScoreMatrix class
    """
    def __init__(self, mat):
        """
        ScoreMatrix constructor
        :param mat: Matrix of substitution scores
        """
        self.mat = mat

    def __str__(self):
        string = ""
        for key in self.mat:
            string += f"{key} {self.mat[key]}\n"
        return string

    def lookup(self, i, j):
        """
            Returns score given by symbol i and j.
        :param i:   Query symbol
        :param j:   Symbol to get score for
        :return:    Entry in sore matrix should it exist, else 0.
        """
        try:
            return self.mat[i][j]
        except KeyError:
            try:
                # Since substitution matrices are symmetric
                return self.mat[j][i]
            except KeyError:
                return 0


class NucleoScoreType(Enum):
    """
    Nucleotide score matrix type enum.
    """
    IDENTITY = 0        # Identity matrix
    NON_UNIFORM = 1     # Tranistion-transversion substitution matrix


class NucleoScoreMatrix(_ScoreMatrix):
    """
    Nucleotide score matrix class.
    """
    def __init__(self, score_type):
        """
        Nucleotide score matrix constructor
        :param score_type: Type of nucleotide scoring matrix to use.
        """
        if score_type in NucleoScoreType:
            self.score_type = score_type
            if score_type is NucleoScoreType.IDENTITY:
                super(NucleoScoreMatrix, self).__init__(_NUCLEO_IDENTITY)
            elif score_type is NucleoScoreType.NON_UNIFORM:
                super(NucleoScoreMatrix, self).__init__(_NUCLEO_NON_UNIFORM)
        else:
            raise TypeError("Param score_type must be of type NucleoScoreType.")


class AminoScoreType(Enum):
    """
    Amino acid score matrix type enum.
    """
    BLOSUM50 = 0
    BLOSUM62 = 1
    BLOSUM80 = 2
    BLOSUM90 = 3
    PAM30 = 4
    PAM70 = 5
    PAM250 = 6


class AminoScoreMatrix(_ScoreMatrix):
    """
    Amino acid score matrix class.
    """
    def __init__(self, score_type):
        """
        Amino acid score matrix constructor.
        :param score_type: Type of amino acid scoring matrix to use.
        """
        if score_type in AminoScoreType:
            self.score_type = score_type
            if score_type is AminoScoreType.BLOSUM62:
                super(AminoScoreMatrix, self).__init__(load(score_type))
        else:
            raise TypeError("Param score_type must be of type AminoScoreType")


_NUCLEO_IDENTITY = {"A": {"A": 1, "T": 0, "C": 0, "G": 0, "U": 0, "-": 0},
                    "T": {"A": 0, "T": 1, "C": 0, "G": 0, "U": 0, "-": 0},
                    "C": {"A": 0, "T": 0, "C": 1, "G": 0, "U": 0, "-": 0},
                    "G": {"A": 0, "T": 0, "C": 0, "G": 1, "U": 0, "-": 0},
                    "U": {"A": 0, "T": 0, "C": 0, "G": 0, "U": 1, "-": 0},
                    "-": {"A": 0, "T": 0, "C": 0, "G": 0, "U": 0, "-": 0}}


_NUCLEO_NON_UNIFORM = {"A": {"A": 1, "T": -1, "C": -0.5, "G": -1, "-": -1},
                       "T": {"A": -1, "T": 1, "C": -1, "G": -0.5, "-": -1},
                       "C": {"A": -0.5, "T": -1, "C": 1, "G": -1, "-": -1},
                       "G": {"A": -1, "T": -0.5, "C": -1, "G": 1, "-": -1},
                       "-": {"A": -1, "T": -1, "C": -1, "G": -1, "-": 1}}

_AMINO_BLOSUM62 = {"A": {"A": 4, "C": 0, "D": -2,"E": -1, "F": -2, "G": 0, "H": -2, "I": -1, "K": -1, "L": -1, "M": -1, "N": -2, "P": -1, "Q": -1, "R": -1, "S": 1, "T": 0, "V": 0, "W": -3, "Y": -2, "*": -4, "-": -1},
                   "C": {"C": 9, "D": -3, "E": -4, "F": -2, "G": -3, "H": -3, "I": -1, "K": -3, "L": -1, "M": -1, "N": -3, "P": -3, "Q": -4, "R": -3, "S": -1, "T": -1, "V": -1, "W": -2, "Y": -2, "*": -4, "-": -1},
                   "D": {"D": 6, "E": 2, "F": -3, "G": -1, "H": -1, "I": -3, "K": -1, "L": -4, "M": -3, "N": 1, "P": -1, "Q": 0, "R": -2, "S": 0, "T": -1, "V": -3, "W": -4, "Y": -3,"*": -4, "-": -1},
                   "E": {"E": 5, "F": -3, "G": -2, "H": 0, "I": -3, "K": 1, "L": -3, "M": -2, "N": 0, "P": -1, "Q": 2, "R": 0, "S": 0, "T": -1, "V": -2, "W": -3, "Y": -2, "*": -4, "-": -1},
                   "F": {"F": 6, "G": -3, "H": -1, "I": 0, "K": -3, "L": 0, "M": 0, "N": -3, "P": -4, "Q": -3, "R": -3, "S": -2, "T": -2, "V": -1, "W": 1, "Y": 3, "*": -4, "-": -1},
                   "G": {"G": 6, "H": -2, "I": -4, "K": -2, "L": -4, "M": -3, "N": 0, "P": -2, "Q": -2, "R": -2, "S": 0, "T": -2, "V": -3, "W": -2, "Y": -3, "*": -4, "-": -1},
                   "H": {"H": 8, "I": -3, "K": -1, "L": -3, "M": -2, "N": 1, "P": -2, "Q": 0, "R": 0, "S": -1, "T": -2, "V": -3, "W": -2, "Y": 2, "*": -4, "-": -1},
                   "I": {"I": 4, "K": -3, "L": 2, "M": 1, "N": -3, "P": -3, "Q": -3, "R": -3, "S": -2, "T": -1, "V": 3, "W": -3, "Y": -1, "*": -4, "-": -1},
                   "K": {"K": 5, "L": -2, "M": -1, "N": 0, "P": -1, "Q": 1, "R": 2, "S": 0, "T": -1, "V": -2, "W": -3, "Y": -2, "*": -4, "-": -1},
                   "L": {"L": 4, "M": 2, "N": -3, "P": -3, "Q": -2, "R": -2, "S": -2, "T": -1, "V": 1, "W": -2, "Y": -1, "*": -4, "-": -1},
                   "M": {"M": 5, "N": -2, "P": -2, "Q": 0, "R": -1, "S": -1, "T": -1, "V": 1, "W": -1, "Y": -1, "*": -4, "-": -1},
                   "N": {"N": 6, "P": -2, "Q": 0, "R": 0, "S": 1, "T": 0, "V": -3, "W": -4, "Y": -2, "*": -4, "-": -1},
                   "P": {"P": 7, "Q": -1, "R": -2, "S": -1, "T": -1, "V": -2, "W": -4, "Y": -3, "*": -4, "-": -1},
                   "Q": {"Q": 5, "R": 1, "S": 0, "T": -3, "V": -2, "W": -2, "Y": -1, "*": -4, "-": -1},
                   "R": {"R": 5, "S": -1, "T": -1, "V": -3, "W": -3, "Y": -2, "*": -4, "-": -1},
                   "S": {"S": 4, "T": 1, "V": -2, "W": -3, "Y": -2, "*": -4, "-": -1},
                   "T": {"T": 5, "V": 0, "W": -2, "Y": -2, "*": -4, "-": -1},
                   "V": {"V": 4, "W": -3, "Y": -1, "*": -4, "-": -1},
                   "W": {"W": 11, "Y": 2, "*": -4, "-": -1},
                   "Y": {"Y": 7, "*": -4, "-": -1},
                   "*": {"*": 1, "-": -1},
                   "-": {"-": -1}}
