from src.seq.Seq import *
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

        # If distribution vector and transition matrix passed
        if "pi" in kwargs and "tm" in kwargs:
            self.tm = kwargs["tm"]
            self.pi = kwargs["pi"]

        # If transition matrix passed but not distribution vector
        elif "tm" in kwargs and "pi" not in kwargs:
            self.tm = kwargs["tm"]
            self.pi = self._rand_pi()

        # If transition matrix not passed but distribution vector is
        elif "tm" not in kwargs and "pi" in kwargs:
            self.tm = self._rand_tm()
            self.pi = kwargs["pi"]
        # Else if neither are passed
        else:
            self.tm = self._rand_tm()
            self.pi = self._rand_pi()

    def _rand_pi(self):
        """
        Randomly generates intitial distriubtion vector.
        :return: Initial distribution vector
        """
        total = 1000    # Total probability
        vec = [0 for _ in self.alpha]   # Probability vector with zeroes

        # Iterate over sequence alphabet
        for i in range(len(self.alpha)):
            # Set proability based on amount left sequence type
            p = 0 if total <= 0 else randint(0, total) // (len(self.alpha) /
                                            (randint(1, 4) if self.seq_type is SeqType.DNA
                                             or self.seq_type is SeqType.RNA
                                             else randint(2, 6)))
            # Last element in distribution vector should equal the
            # equal the total amount left
            # Else whatever p was calculated to be
            # Both divided by 100-
            vec[i] = (total if i == len(self.alpha) -1 else p) / 1000
            # Subtract from total probability
            total -= p
        # Return probability distribution vector
        return vec

    def _rand_tm(self):
        """
        Generates random transformation matrix.
        :return: Transition matrix.
        """
        # Generate transition matrix with all zeroes
        mat = [[0 for __ in self.alpha] for _ in self.alpha]

        # Iterate over length of alphabet
        for i in range(0, len(self.alpha)):
            # Initialise total probability to 10000
            total = 10000
            # Iterate over length of alphabet
            for j in range(0, len(self.alpha)):
                # Current probability p is zero if there is none left
                # Else it is random integer between 0 and total
                # Floor-divided by a value dictated by sequence type
                #   and hence sequence alphabet type
                p = 0 if total <= 0 else randint(0, total)\
                                         // (len(self.alpha) /
                                             (randint(1,4) if self.seq_type is SeqType.DNA
                                              or self.seq_type is SeqType.RNA
                                              else randint(2,6)))
                # If we are at the last letter of the alphabet
                # Then the probability element
                # should be all of what's left / 10000

                # Else it is the value of p / 10000
                mat[i][j] = (total if j == len(self.alpha)-1 else p) / 10000
                # Subtract the value of p from the total
                total -= p
        # Return the transition matrix
        return mat

    def sequence(self, length=-1):
        """
        Generates biosequence
        :param length: Length of sequence to generate
        :return:       Generated sequence
        """
        seq = "" # Sequence to be generated
        index = 0 # Current index
        total_prob = 0 # Total probability initialised to zero
        p = random() # Random value
        seq_length = length # Sequence length

        # If length not specified generate sequence between 1 and 30 residues
        if length == -1:
            seq_length = randint(1, 30)
        # If length specified is incorrect type do the same
        elif not isinstance(length, int):
            seq_length = randint(1, 30)

        # Iterate over probability vector
        for i in range(len(self.pi)):
            # Sum probabilities
            total_prob += self.pi[i]
            # If probability greater than random
            if total_prob > p:
                # Set index
                index = i

        # Loop while residues left
        while len(seq) < seq_length:
            # Get random probability
            p = random()
            # Set sum of probabilities to zero
            total_prob = 0
            # Append residue using seq alphabet and index
            seq += self.alpha[index]
            # Iterate over rows of transition matrix
            for j in range(len(self.tm[index])):
                # Sum probabilities
                total_prob += self.tm[index][j]
                # If sum greater than random
                if total_prob > p:
                    # Transition to j
                    index = j
                    # End loop
                    break

        # Return generated sequence as a Seq object
        return Seq(seq, self.seq_type)



