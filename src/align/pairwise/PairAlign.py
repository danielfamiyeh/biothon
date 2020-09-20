from src.seq.Seq import Seq

class PairAlign:
    """
    Pairwise alignment class for printing pairwise alignments
        to terminal.
    """
    def __init__(self, seq1, seq2, pipes, score):
        self.seq1 = seq1.split("\n")
        self.seq2 = seq2.split("\n")
        self.pipes = "".join(pipes).split("\n")
        self.shift_length = len(seq1.id) if len(seq1.id) >\
                                              len(seq2.id) else len(seq2.id)
        self.score = score

        mismatches = 0
        no_gap = 0
        # Iterate over both sequences
        for char1, char2 in zip(self.seq1, self.seq2):
            if char1 != "-" and char2 != "-":
                # Increment mismatch count
                mismatches += int(char1 != char2)
                # Increment not gapped count
                no_gap += 1

        # Calculate distance to 3 dp
        self.dist = round((mismatches / (no_gap if no_gap > 0 else 1) * 100), 3)

    def __str__(self):
        return ''.join([seq1.id + (' ' * (self.shift_length - len(seq1.id))) +
                        ("\t\t" if self.shift_length > 0 else "") + seq1.seq + "\n"
                        + (' ' * self.shift_length) +
                        ("\t\t" if self.shift_length > 0 else "") + pipes + "\n"
                        + seq2.id + (' ' * (self.shift_length - len(seq2.id))) +
                        ("\t\t" if self.shift_length > 0 else "") + seq2.seq + "\n"
                        for seq1, seq2, pipes in zip(self.seq1, self.seq2, self.pipes)])

    def __repr__(self):
        return f"Alignment Score: {self.score}\n" + self.__str__()