from src.seq.Seq import Seq

class PairAlign:
    """
    Pairwise alignment class for printing pairwise alignments
        to terminal.
    """
    def __init__(self, seq1, seq2, pipes):
        self.seq1 = seq1.split("\n")
        self.seq2 = seq2.split("\n")
        self.pipes = "".join(pipes).split("\n")
        self.shift_length = len(seq1.name) if len(seq1.name) >\
                                              len(seq2.name) else len(seq2.name)

    def __str__(self):
        return ''.join([seq1.name + (' ' * (self.shift_length - len(seq1.name))) + "\t\t" + seq1.seq + "\n"
                        + (' ' * self.shift_length) + "\t\t" + pipes + "\n"
                        + seq2.name + (' ' * (self.shift_length - len(seq2.name))) + "\t\t" + seq2.seq + "\n\n"
                        for seq1, seq2, pipes in zip(self.seq1, self.seq2, self.pipes)])