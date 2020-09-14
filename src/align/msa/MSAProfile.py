from functools import reduce


class MSAProfile:
    def __init__(self, seq, score_matrix):
        """
        Multiple sequence alignment profile constructor
        :param seq: Initial sequence
        :param score_matrix: Score matrix for subsequent alignment
        """
        
        # Sequence type assigment
        if seq.seq_type == SeqType.DNA:
            self.alpha = alphabet_dna
        elif seq.seq_type == SeqType.RNA:
            self.alpha = alphabet_rna
        elif seq.seq_type == SeqType.PROTEIN:
            self.alpha = alphabet_protein_1
        else:
            raise TypeError("Param seq must be have valid SeqType")

        # Required (and updated) in order to align correctly
        self.len_longest_name = len(seq.name)

        self.profile = {c: [] for c in self.alpha}
        self.seq_list = [seq]
        self.consensus = seq
        self.score_matrix = score_matrix

        # Initialising profile
        for k, v in self.profile.items():
            for c in seq.seq:
                v.append(1 if c == k else 0)

    def get_alignments(self):
        """
        Returns multiple sequence alignment as a string.
        :return: String containing alignment data in clustal format
        """
        tab = '\t'
        string = ""
        char_index = 0
        base_counts = [0 for _ in self.seq_list] # Rolling alignment count

        # Generates alignment string with max 60 residues per line
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
        '''
        Returns profile as a string
        :return: Character profile string of multiple sequence alignment
        '''
        return "\n".join([f"{k} {v}" for k, v in self.profile.items()])

    def align(self, **kwargs):
        '''
        ClustalW multiple sequence alignment algorithm
        :param kwargs:
            seq (Seq): Sequence to add to alignment
            profile (MSAProfile): Multiple sequences to add to alignment
        :return: None
        '''
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

            # Matrix to keep track of predecessors in alignment
            predecessors = [[0 for _ in range(profile_len)]
                            for __ in range(query_len)]

            # Populates set of letters in MSA
            for i in range(1, len(matrix)):
                for j in range(1, len(matrix[0])):
                    letters = set()
                    for c in self.alpha:
                        if self.profile[c][j - 1] != 0:
                            letters.add(c)

                    try:
                        # TODO: match weight with letter - This doesn't feel like it should work.
                        weight = reduce(lambda x, y: x * y, [s.weight for s in self.seq_list if s[j - 1] in letters])
                    except TypeError:
                        weight = 1

                    if is_seq:  # Then we are performing Seq-MultiSeq alignment
                        
                        # TODO: replace with for/for-in to solve weighting issue
                        current = sum(self.profile[self.alpha[k]][j - 1] *
                                      self.score_matrix.lookup(query[i - 1], self.alpha[k])
                                      for k in range(len(self.alpha))) * weight

                    else:   # Then we are performing MultiSeq-MultiSeq alignment
                        current = 0

                        # Iterate over all sequences in query MultiSeq
                        for query_seq in query.seq_list:
                            current += sum(self.profile[self.alpha[k]][j - 1] *
                                           self.score_matrix.lookup(query_seq[i - 1], self.alpha[k])
                                           for k in range(len(self.alpha)))
                        current /= query_len
                        current *= weight

                    scores = [matrix[i - 1][j - 1] + current,
                              matrix[i][j - 1] + -1,
                              matrix[i - 1][j] + -1]

                    # Assign value and predecessor
                    argmax = max(scores)
                    matrix[i][j] = argmax
                    predecessors[i - 1][j - 1] = scores.index(argmax)

            # Set i, j values to end of matrix
            i = len(predecessors)
            j = len(predecessors[0])

            # Aligned sequences, query sequence and list of query sequences.
            aligned_seqs = ["" for seq in self.seq_list]
            aligned_query = ""
            aligned_queries = None if is_seq else ["" for _ in range(query_len)]

            # Backtracking
            while i > 0 or j > 0:
                if i == 0:  # We have reached the first character of original MSA
                    j -= 1
                    # Iterate from length of query sequence backwards to zero
                    for count in range(j, -1, -1):
                        j -= 1
                        if aligned_queries is None:
                            # Then we are performing Seq-MultiSeq alignment
                            aligned_query += "-"    # So add gaps to aligned query seq
                        else:
                            # Then we are perfoming MultiSeq-MultiSeq alignment
                            for k in range(query_len):
                                # So iterate over all query sequences
                                aligned_queries[k] += "-"   # And add gaps to them

                        for k in range(len(self.seq_list)):
                            # In either case add the letters from the MSA to new
                            #   aligned_seqs list
                            aligned_seqs[k] += self.seq_list[k][count]

                elif j == 0:    # Then we are at the beginning of the query sequence
                    i -= 1
                    for count in range(i, -1, -1):
                        i -= 1
                        if aligned_queries is None:
                            # Then we are performing Seq-MultiSeq alignment
                            aligned_query += query[count]
                        else:
                            # Then we are performing MultiSeq-MultiSeq alignment
                            for k in range(query_len):
                                # Add rest of letters to aligned_queries list
                                aligned_queries[k] += query[k][count]

                        for k in range(len(self.seq_list)):
                            # Add gaps to aligned db sequences
                            aligned_seqs[k] += "-"

                # Diagonal case
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

                # Left (query gap) case
                elif predecessors[i - 1][j - 1] == 1:
                    j -= 1
                    if aligned_queries is None:
                        aligned_query += "-"
                    else:
                        for k in range(query_len):
                            aligned_queries[k] += "-"
                    for k in range(len(self.seq_list)):
                        aligned_seqs[k] += self.seq_list[k][j]

                # Above (db seq gap) case
                else:
                    i -= 1
                    if aligned_queries is None:
                        aligned_query += query[i]
                    else:
                        for k in range(query_len):
                            aligned_queries[k] += query.seq_list[k][i]
                    for k in range(len(self.seq_list)):
                        aligned_seqs[k] += "-"

            # Reverse all aligned sequences
            for k in range(len(self.seq_list)):
                self.seq_list[k].seq = aligned_seqs[k][::-1]

            # If this is Seq-MultiSeq alignment
            if is_seq:
                # Reverse query string and add new Seq to seq_list
                self.seq_list.append(Seq(aligned_query[::-1], query.seq_type,
                                         name=query.name, label=query.label, weight=query.weight))
            else:
                # If this is MultiSeq-MultiSeq alignment
                for i, string in enumerate(aligned_queries):
                    # Add and reverse all new sequences as Seq objects to seq_list
                    self.seq_list.append(Seq(string[::-1], query.seq_list[i].seq_type, name=query.seq_list[i].name,
                                             label=query.seq_list[i].label, weight=query.seq_list[i].weight))

        else:
            # Raise error if alignment args are is invalid
            raise TypeError("Can only perform sequence-profile and profile-profile alignments.")

        self.reprofile()    # Re-calculate profile distribution

    def reprofile(self):
        """
        Recalculates character distributions in profile.
        :return: None
        """
        self.profile = {c: [0 for _ in range(len(self.seq_list[0]))]
                        for c in self.alpha} # New empty profile

        # Character counter for every column in profile
        column_count = [{k: 0 for k in self.alpha} for _ in range(len(self.seq_list[0]))]

        # Gap counter for every column in profile
        column_gap = [0 for _ in range(len(self.seq_list[0]))]

        # Consensus tuples for every column in profile
        consensuses = [("", 0) for _ in range(len(self.seq_list[0]))]

        for seq in self.seq_list:   # Iterate over all sequences
            for i, c in enumerate(seq):
                # Enumerate every letter in sequence
                if c == "-":
                    # Increment gap counter if character is a gap
                    column_gap[i] += 1
                else:
                    # Increment character counter if character is not a gap
                    column_count[i][c] += 1

            # Update longest sequence name if necessary
            if len(seq.name) > self.len_longest_name:
                self.len_longest_name = len(seq.name)

        # Enumerate sequence alphabet
        for j, c in enumerate(self.alpha):
            # Iterate over seq_list
            for i in range(len(self.seq_list[0])):
                # Calculate divisor for column's distribution
                divisor = len(self.seq_list) - column_gap[i]
                # Divide through every entry by their corresponding divisor
                self.profile[c][i] = column_count[i][c] / (divisor if divisor > 0 else 1)

        # Iterate over key/values in profile
        for k, v in self.profile.items():
            # Round entries to two decimal places
            self.profile[k] = list(map(lambda x: round(x, 2), v))

        # Iterate over every character in sequence alphabet
        for char in self.alpha:
            # Enumerate over character distributions (by row)
            for i, val in enumerate(self.profile[char]):
                if val > consensuses[i][1]:
                    # Update consensus character if a new one is found
                    consensuses[i] = (char, val)

        # Update consensus string
        self.consensus = "".join([c[0] for c in consensuses])

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"Biothon ClustalW Multiple Sequence Alignment\n\n{self.get_alignments()}"

    def __len__(self):
        return len(self.seq_list)
