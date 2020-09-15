class PairAlignment:
    def __init__(self, *args):
        """
        Pairwise alignment object constructor.
        :param args:    Two aligned strings and their alignment score
        :param kwargs:  s1_name, s2_name representing the names
                            given to the sequences during construction
        """

        mismatches, no_gap = 0, 0

        # Assign sequence data
        self.seq1, self.seq2, self.score = args

        self.seq_type = self.seq1.seq_type

        # Length of longest sequence name
        self.len_longest_id = max(len(self.seq1.id), len(self.seq2.id))

        # Iterate over both sequences
        for char1, char2 in zip(self.seq1, self.seq2):
            if char1 != "-" and char2 != "-":
                # Increment mismatch count
                mismatches += int(char1 != char2)
                # Increment not gapped count
                no_gap += 1

        # Calculate distance to 3 dp
        self.dist = round((mismatches/(no_gap if no_gap > 0 else 1) * 100), 3)

        self.string = self._form_string()

    def _form_string(self):
        """
        Forms a string to represent the alignment using pipes between like-bases.
        :return: Aligned string
        """
        string = ""
        tab = '\t' if self.len_longest_id > 0 else ''

        # Prefixes i.e. Seq IDs and tabbed spacing to ensure aligned seqs match up
        before_seq1 = self.seq1.id + ' ' * (self.len_longest_id - len(self.seq1.id)) + tab
        before_pipes = ' ' * self.len_longest_id + tab
        before_seq2 = self.seq2.id + ' ' * (self.len_longest_id - len(self.seq2.id)) + tab

        splits = [[], [], []]
        prefixes = [before_seq1, before_pipes, before_seq2]
        for i in range(3):
            count = - 1
            for j in range(len(self.seq1)):
                if j % 60 == 0:
                    splits[i].append([prefixes[i]])
                    count += 1
                splits[i][count].append(
                    self.seq1[j] if i == 0 else
                    ("|" if self.seq1[j] == self.seq2[j] else " ") if i == 1
                    else self.seq2[j]

                )

        for i in range(len(splits[0])):
            string += ''.join(splits[0][i]) + "\n"
            string += ''.join(splits[1][i]) + "\n"
            string += ''.join(splits[2][i]) + "\n\n"
        return string

    def __str__(self):
        return self.string

    def __repr__(self):
        return f"Pairwise Alignment | Score: {self.score} | " \
               f"Distance: {self.dist}%\n\n{self.string}"

    def __len__(self):
        return len(self.seq1)
