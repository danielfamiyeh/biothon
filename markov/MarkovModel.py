from seq.Seq import *
from random import random, randint


class MarkovModel:
    """
    MarkovModel class for probabilistic sequence generation.
    """

    def __init__(self, **kwargs):
        """
        MarkovModel constructor
        :param kwargs: type: Sequence type
                       pi:   Initial distribution vector
                       tm:   Transition matrix
        """
        self.seq_type = kwargs.get("type", SeqType.DNA)
        self.alpha = alphabet_dna if self.seq_type is SeqType.DNA else\
                     alphabet_rna if self.seq_type is SeqType.RNA else\
                    alphabet_protein_1

        if "pi" in kwargs and "tm" in kwargs:
            self.tm = kwargs["tm"]
            self.pi = kwargs["pi"]

        elif "tm" in kwargs and "pi" not in kwargs:
            self.tm = kwargs["tm"]
            self.pi = self._rand_pi()

        elif "tm" not in kwargs and "pi" in kwargs:
            self.tm = self._rand_tm()
            self.pi = kwargs["pi"]
        else:
            self.tm = self._rand_tm()
            self.pi = self._rand_pi()

    def _rand_pi(self):
        """
        Randomly generates intitial distriubtion vector.
        :return: Initial distribution vector
        """
        total = 1000
        vec = [0 for _ in self.alpha]
        for i in range(len(self.alpha)):
            p = 0 if total <= 0 else randint(0, total) // (len(self.alpha) /
                                            (randint(1, 4) if self.seq_type is SeqType.DNA
                                             or self.seq_type is SeqType.RNA
                                             else randint(2, 6)))
            vec[i] = (total if i == len(self.alpha) -1 else p) / 1000
            total -= p
        return vec

    def _rand_tm(self):
        """
        Generates random transformation matrix.
        :return: Transition matrix.
        """
        mat = [[0 for __ in self.alpha] for _ in self.alpha]
        for i in range(0, len(self.alpha)):
            total = 10000
            for j in range(0, len(self.alpha)):
                p = 0 if total <= 0 else randint(0, total)\
                                         // (len(self.alpha) /
                                             (randint(1,4) if self.seq_type is SeqType.DNA
                                              or self.seq_type is SeqType.RNA
                                              else randint(2,6)))
                mat[i][j] = (total if j == len(self.alpha)-1 else p) / 10000
                total -= p
        return mat

    def sequence(self, length=-1):
        """
        Generates biosequence
        :param length: Length of sequence to generate
        :return:       Generated sequence
        """
        seq = ""
        index = 0
        total_prob = 0
        p = random()
        seq_length = length

        if length == -1:
            seq_length = randint(1, 30)
        elif not isinstance(length, int):
            seq_length = randint(1, 30)

        for i in range(len(self.pi)):
            total_prob += self.pi[i]
            if total_prob > p:
                index = i

        while len(seq) < seq_length:
            p = random()
            total_prob = 0
            seq += self.alpha[index]
            for j in range(len(self.tm[index])):
                total_prob += self.tm[index][j]
                if total_prob > p:
                    index = j
                    break

        return Seq(seq, self.seq_type)



