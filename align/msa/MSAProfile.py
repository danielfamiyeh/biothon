from functools import reduce
from tree.GuideTree import *


class MSAProfile:
    def __init__(self, seq, score_matrix):
        """
        Multiple sequence alignment profile constructor
        :param seq: Initial sequence
        :param score_matrix: Score matrix for subsequent alignment
        """
        if seq.seq_type == SeqType.DNA:
            self.alpha = alphabet_dna
        elif seq.seq_type == SeqType.RNA:
            self.alpha = alphabet_rna
        elif seq.seq_type == SeqType.PROTEIN:
            self.alpha = alphabet_protein_1
        else:
            raise TypeError("Param seq must be have valid SeqType")

        self.len_longest_name = len(seq.name)
        self.profile = {c: [] for c in self.alpha}
        self.seq_list = [seq]
        self.consensus = seq
        self.score_matrix = score_matrix

        for k, v in self.profile.items():
            for c in seq.seq:
                v.append(1 if c == k else 0)

    def get_alignments(self):
        """
        :return: String containing alignment data in clustal format
        """
        tab = '\t'
        string = ""
        char_index = 0
        base_counts = [0 for _ in self.seq_list]
        while char_index < len(self.seq_list[0]):
            for i, s in enumerate(self.seq_list):
                if char_index <= len(s):
                    string += f"{s.name}{' '*(self.len_longest_name-len(s.name)) + tab*2}"
                    for j in range(char_index, char_index+60):
                        if j >= len(s):
                            break
                        else:
                            if s[j] != "-":
                                base_counts[i] += 1
                            string += s[j]
                    string += f" {base_counts[i]}\n"
            char_index += 60
            string += "\n"
        return string

    def get_profile(self):
        return "\n".join([f"{k} {v}" for k, v in self.profile.items()])

    def align(self, **kwargs):
        # Sequence-Profile Alignment
        if "seq" in kwargs or "profile" in kwargs:
            is_seq = kwargs.get("seq", False)
            query = kwargs["seq" if is_seq else "profile"]

            query_len = len(query)
            profile_len = len(self.profile[self.alpha[0]])

            # Generating alignment matrix
            matrix = [[0 if i == 0 and j == 0 else 0 if i > 0 and j > 0 else
                      -1 * j if i == 0 and j >= 0 else -1*i
                      for j in range(profile_len + 1)]
                      for i in range(query_len + 1)]

            predecessors = [[0 for _ in range(profile_len)]
                            for __ in range(query_len)]

            for i in range(1, len(matrix)):
                for j in range(1, len(matrix[0])):
                    letters = set()
                    for c in self.alpha:
                        if self.profile[c][j - 1] != 0:
                            letters.add(c)

                    try:
                        weight = reduce(lambda x, y: x * y, [s.weight for s in self.seq_list if s[j - 1] in letters])
                    except TypeError:
                        weight = 1

                    if is_seq:
                        current = sum(self.profile[self.alpha[k]][j - 1] *
                                      self.score_matrix.lookup(query[i - 1], self.alpha[k])
                                      for k in range(len(self.alpha))) * weight
                    else:
                        current = 0
                        for q in query.seq_list:
                            current += sum(self.profile[self.alpha[k]][j - 1] *
                                           self.score_matrix.lookup(q[i - 1], self.alpha[k])
                                           for k in range(len(self.alpha)))
                        current /= query_len
                        current *= weight
                    scores = [matrix[i - 1][j - 1] + current, matrix[i][j - 1] + -1, matrix[i - 1][j] + -1]
                    argmax = max(scores)
                    matrix[i][j] = argmax
                    predecessors[i - 1][j - 1] = scores.index(argmax)

            i = len(predecessors)
            j = len(predecessors[0])

            aligned_seqs = ["" for seq in self.seq_list]
            aligned_query = ""
            aligned_queries = None if is_seq else ["" for _ in range(query_len)]

            # Backtracking
            while i > 0 or j > 0:
                if i == 0:
                    j -= 1
                    for count in range(j, -1, -1):
                        j -= 1
                        if aligned_queries is None:
                            aligned_query += "-"
                        else:
                            for k in range(query_len):
                                aligned_queries[k] += "-"
                        for k in range(len(self.seq_list)):
                            aligned_seqs[k] += self.seq_list[k][count]

                elif j == 0:
                    i -= 1
                    for count in range(i, -1, -1):
                        i -= 1
                        if aligned_queries is None:
                            aligned_query += query[count]
                        else:
                            for k in range(query_len):
                                aligned_queries[k] += query[k][count]
                        for k in range(len(self.seq_list)):
                            aligned_seqs[k] += "-"

                # diagonal
                elif predecessors[i - 1][j - 1] == 0:
                    i -= 1
                    j -= 1
                    if aligned_queries is None:
                        aligned_query += query[i]
                    else:
                        for k in range(query_len):
                            aligned_queries[k] += query.seq_list[k][i]
                    for k in range(len(self.seq_list)):
                        aligned_seqs[k] += self.seq_list[k][j]

                # left
                elif predecessors[i - 1][j - 1] == 1:
                    j -= 1
                    if aligned_queries is None:
                        aligned_query += "-"
                    else:
                        for k in range(query_len):
                            aligned_queries[k] += "-"
                    for k in range(len(self.seq_list)):
                        aligned_seqs[k] += self.seq_list[k][j]

                # up
                else:
                    i -= 1
                    if aligned_queries is None:
                        aligned_query += query[i]
                    else:
                        for k in range(query_len):
                            aligned_queries[k] += query.seq_list[k][i]
                    for k in range(len(self.seq_list)):
                        aligned_seqs[k] += "-"

            for k in range(len(self.seq_list)):
                self.seq_list[k].seq = aligned_seqs[k][::-1]

            if is_seq:
                self.seq_list.append(Seq(aligned_query[::-1], query.seq_type,
                                         name=query.name, label=query.label, weight=query.weight))
            else:
                for i, string in enumerate(aligned_queries):
                    self.seq_list.append(Seq(string[::-1], query.seq_list[i].seq_type, name=query.seq_list[i].name,
                                             label=query.seq_list[i].label, weight=query.seq_list[i].weight))

        else:
            raise TypeError("Can only perform sequence-profile and profile-profile alignments.")

        self.reprofile()

    def reprofile(self):
        """
        Recalculats character distributions in profile.
        :return: None
        """
        self.profile = {c: [0 for _ in range(len(self.seq_list[0]))]
                        for c in self.alpha}

        column_count = [{k: 0 for k in self.alpha} for _ in range(len(self.seq_list[0]))]
        column_gap = [0 for _ in range(len(self.seq_list[0]))]
        consensuses = [("", 0) for _ in range(len(self.seq_list[0]))]

        for seq in self.seq_list:
            for i, c in enumerate(seq):
                if c == "-":
                    column_gap[i] += 1
                else:
                    column_count[i][c] += 1
            if len(seq.name) > self.len_longest_name:
                self.len_longest_name = len(seq.name)

        for j, c in enumerate(self.alpha):
            for i in range(len(self.seq_list[0])):
                divisor = len(self.seq_list) - column_gap[i]
                self.profile[c][i] = column_count[i][c] / (divisor if divisor > 0 else 1)

        for k, v in self.profile.items():
            self.profile[k] = list(map(lambda x: round(x, 2), v))

        for char in self.alpha:
            for i, val in enumerate(self.profile[char]):
                if val > consensuses[i][1]:
                    consensuses[i] = (char, val)
        self.consensus = "".join([c[0] for c in consensuses])

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"Biothon ClustalW Multiple Sequence Alignment\n\n{self.get_alignments()}"

    def __len__(self):
        return len(self.seq_list)
