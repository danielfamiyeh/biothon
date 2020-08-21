from math import ceil
from align.pair.PairAlignment import PairAlignment
from seq.Seq import *
from collections import namedtuple
from graph.Graph import Graph

class Pair:
    def __init__(self, val):
        """
        Pair object constructor, used to represent entry in matrix
        :param val: Value of entry in matrix.
        """
        self.val = val
        self.pred = None    # Predecessor 'pointer' makes backtracking easier.

    def __repr__(self):
        return str(self.val)

    def __sub__(self, other):
        if isinstance(other, int):
            return self.val - other


class PairAligner:
    def __init__(self, score_matrix, match=1, mismatch=0, gap_open=1):
        """
        Pairwise aligner object constructor.
        :param score_matrix:    Scoring matrix for alignments.
        :param match:           Score for a match between sequences.
        :param mismatch:        Mismatch penalty
        :param gap_open:        Gap penalty
        """
        self.score_matrix = score_matrix
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open

    def needle(self, s1, s2):
        """
        Standard Needleman-Wunsch algorithm for global pairwise alignment.
        :param s1:              First query sequence
        :param s2:              Second query sequence
        :return:                Returns 3-tuple of aligned sequences and score.
        """
        if isinstance(s1, Seq) and isinstance(s2, Seq):
            score = 0

            rows = len(s1) + 1
            cols = len(s2) + 1
            matrix = []

            # Generate scoring matrix
            for i in range(rows):
                matrix.append([])
                for j in range(cols):
                    matrix[i].append(Pair(0))
                    if i == 0:
                        matrix[i][j].val = - self.gap_open * j
                    elif j == 0:
                        matrix[i][j].val = - self.gap_open * i

            # Traverse matrix and update scores
            for i in range(1, rows):
                for j in range(1, cols):
                    scores = [(matrix[i][j-1]-self.gap_open),
                              (matrix[i-1][j]-self.gap_open)]

                    tentative_match = matrix[i - 1][j - 1].val
                    tentative_match += self.match if s1[i - 1] == s2[j - 1] else -self.mismatch

                    scores.insert(0, tentative_match)

                    matrix[i][j].val = max(scores)
                    matrix[i][j].pred = scores.index(matrix[i][j].val)

            i = rows - 1
            j = cols - 1

            aligned_s1 = ""
            aligned_s2 = ""

            k = 0

            while i > 0 or j > 0:
                if i == 0:
                    j -= 1
                    for count in range(j, -1, -1):
                        k += 1
                        j -= 1
                        aligned_s1 += "-"
                        aligned_s2 += s2[count]
                        score -= self.gap_open

                elif j == 0:
                    i -= 1
                    for count in range(i, -1, -1):
                        k += 1
                        i -= 1
                        aligned_s1 += s1[count]
                        aligned_s2 += "-"
                        score -= self.gap_open

                elif matrix[i][j].pred == 0:
                    i -= 1
                    j -= 1
                    k = 0

                    aligned_s1 += s1[i]
                    aligned_s2 += s2[j]

                    score += self.score_matrix.lookup(s1[i], s2[j])

                elif matrix[i][j].pred == 1:
                    j -= 1
                    k += 1

                    aligned_s1 += "-"
                    aligned_s2 += s2[j]

                    score -= self.gap_open

                else:
                    i -= 1
                    k += 1

                    aligned_s1 += s1[i]
                    aligned_s2 += "-"

                    score -= self.gap_open

            return PairAlignment(aligned_s1[::-1], aligned_s2[::-1], score,
                                 s1_name=s1.name, s2_name=s2.name)
        else:
            raise TypeError("Params s1 and s2 must be of type Seq.")

    def smith(self, s1, s2):
        """
        Standard Smith-Waterman algorithm for pairwise local alignment
        :param s1:              First query sequence
        :param s2:              Second query sequence
        :return:                Returns 3-tuple of aligned sequences and score.
        """

        if isinstance(s1, Seq) and isinstance(s2, Seq):
            score = 0

            rows = len(s1) + 1
            cols = len(s2) + 1

            matrix = []

            # Generate zero matrix
            for i in range(rows):
                matrix.append([Pair(0) for _ in range(cols)])

            for i in range(1, rows):
                for j in range(1, cols):
                    scores = [(matrix[i][j-1] - self.gap_open), (matrix[i-1][j] - self.gap_open), 0]

                    tentative_match = matrix[i-1][j-1].val
                    tentative_match += self.match if s1[i-1] == s2[j-1] else -self.mismatch

                    scores.insert(0, tentative_match)

                    matrix[i][j].val = max(scores)
                    matrix[i][j].pred = scores.index(matrix[i][j].val)

            highest = 0
            max_row = 0
            max_col = 0

            s1_local = ""
            s2_local = ""

            k = 0

            for i, row in enumerate(matrix):
                for j, col in enumerate(row):
                    if highest < col.val:
                        highest = col.val
                        max_row = i
                        max_col = j

            while matrix[max_row][max_col].val > 0:
                if matrix[max_row][max_col].pred == 0:
                    max_row -= 1
                    max_col -= 1
                    k = 0

                    s1_local += s1[max_row]
                    s2_local += s2[max_col]

                    score += self.score_matrix.lookup(s1[max_row], s2[max_col])

                elif matrix[max_row][max_col].pred == 1:
                    max_col -= 1
                    k += 1

                    s1_local += "-"
                    s2_local += s2[max_col]

                    score -= self.gap_open

                else:
                    max_row -= 1
                    k += 1

                    s1_local += s1[max_row]
                    s2_local += "-"

                    score -= self.gap_open

            return PairAlignment(s1_local[::-1], s2_local[::-1], score,
                                 s1_name=s1.name, s2_name=s2.name)
        else:
            raise TypeError("Params s1 and s2 must be of type Seq.")

    def __needle_score(self, s1, s2):
        """
        Optimised Needleman-Wunsch that computes scores only.
        Score is computed in O(min{|s1|, |s2|}) (linear) space.
        :param s1:              First query sequence
        :param s2:              Second query sequence
        :return:                Returns last line of Needleman-Wunsch alignment.
        """
        rows = len(s1) + 1
        cols = len(s2) + 1

        matrix = [[0 if i == 1 else j * -self.gap_open for j in range(cols)] for i in range(2)]

        for i in range(1, rows):
            matrix[1][0] = i * -self.gap_open
            for j in range(1, cols):
                scores = [(matrix[1][j - 1] - self.gap_open), (matrix[0][j] - self.gap_open)]
                tentative_match = matrix[0][j-1]
                tentative_match += self.match if s1[i - 1] == s2[j - 1] else -self.mismatch

                scores.insert(0, tentative_match)
                matrix[1][j] = max(scores)

            matrix[0] = matrix[1].copy()
            matrix[1] = [0 for _ in range(cols)]

        return matrix[0]

    def _hirschberg_global(self, s1, s2, score):
        """
        Internal hirschberg algorithm that gets called recursively.
        :param s1: Sequence to align
        :param s2: Sequence to align
        :param score: Running total of score
        :return: 3-tuple containing two aligned strings and a score
        """
        str1 = ""
        str2 = ""

        if len(s1) == 0:
            for char in s2:
                str1 += "-"
                str2 += char
                score -= self.gap_open

        elif len(s2) == 0:
            for char in s1:
                str1 += char
                str2 += "-"
                score -= self.gap_open

        elif len(s1) == 1 or len(s2) == 1:
            t = self.needle(s1, s2)
            str1 += t.seq1
            str2 += t.seq2
            score += t.score

        else:
            y_mid = ceil(len(s1) // 2)

            top = self.__needle_score(s1[0:y_mid], s2)
            bot = self.__needle_score((s1[y_mid:])[::-1], s2[::-1])[::-1]
            sum_vec = [top[i] + bot[i] for i in range(len(s2))]
            x_mid = sum_vec.index(max(sum_vec))

            t1 = self._hirschberg_global(Seq(s1[0:y_mid], s1.seq_type), Seq(s2[0: x_mid], s2.seq_type), score)
            t2 = self._hirschberg_global(Seq(s1[y_mid:], s1.seq_type), Seq(s2[x_mid:], s2.seq_type), score)

            str1 += t1[0] + t2[0]
            str2 += t1[1] + t2[1]
            score += t1[2] + t2[2]

        return str1, str2, score

    def hirschberg_global(self, s1, s2):
        if isinstance(s1, Seq) and isinstance(s2, Seq):
            score = 0
            return PairAlignment(*self._hirschberg_global(s1, s2, score))
        else:
            raise TypeError("Param s1 and s2 must be of type Seq.")

    def blast(self, **kwargs):
        """
        BLAST1 local alignment algorithm.
        :param kwargs: Query string, database string list, k-length, threshold, cutoff score and
                        minimum score for a HSP to be considered statistically significant.
        :return: Sorted list of HSP namedtuples.
        """

        query = kwargs["query"]
        db = kwargs["db"]
        ktup = kwargs["ktup"]
        threshold = kwargs["threshold"]
        cutoff = kwargs["cutoff"]
        min_score = kwargs["min_score"]

        highest_score = 0
        neighbourhood = {}
        hits = {i: set() for i in range(len(db))}
        hsp = []

        alpha = alphabet_dna if query.seq_type is SeqType.DNA else alphabet_rna \
            if query.seq_type is SeqType.RNA else alphabet_protein_1

        # Populating k-word neighbourhood
        for i in range(len(query) - ktup):
            for c in alpha:
                kword = c + query[i+1:ktup]
                query_sub = query[i:ktup]

                match = 0
                for j in range(len(query_sub)):
                    match += self.score_matrix.lookup(query_sub[j], kword[j])

                if match > threshold:
                    if i not in neighbourhood:
                        neighbourhood[i] = []
                    neighbourhood[i].append(kword)

        # Getting seed hits
        for i, seq in enumerate(db):
            for key, value in neighbourhood.items():
                for j, neighbour in enumerate(value):
                    k = 0
                    while True:
                        if k > len(seq):
                            break
                        match = seq[k: k+ktup].find(neighbour)
                        if match != -1:
                            hits[i].add((key, j, k + match))
                        k += 1

        HSPTuple = namedtuple("HSP", [f"seq_label", "subseq", "score", "query_index", "seq_index"])

        # Getting high scoring pairs
        for i, hit_list in hits.items():
            if len(hit_list) > 0:
                for j, seed_hit in enumerate(hit_list):
                    left = 0
                    right = 0

                    while True:
                        query_sub = query[seed_hit[0] - left: seed_hit[0] + ktup + right]
                        seq_sub = db[i][seed_hit[2] - left: seed_hit[2] + ktup + right]

                        match = 0
                        for k in range(min(len(seq_sub), len(query_sub))):
                            match += self.score_matrix.lookup(query_sub[k], seq_sub[k])

                        if match >= min_score and match > highest_score:
                            hsp.append(HSPTuple((db[i].label if db[i].label != "" else i),
                                                seq_sub, match, seed_hit[0] - left,
                                                seed_hit[2] - left))

                            highest_score = match

                        elif match < highest_score - cutoff:
                            break

                        elif min(seed_hit[0] - left, seed_hit[2] - left) < 0 and\
                            (seed_hit[2] + ktup + right > len(db[i]) or
                             seed_hit[0] + ktup + right > len(query)):
                            break

                        left += 1
                        right += 1

        hsp = sorted(hsp, key=lambda x: x.score, reverse=True)
        return hsp

    def _banded_sw(self, runs, query, db, w):
        """
        Banded Smith-Waterman algorithm. Performs local alignment that is
        constrained to a certain distance from a diagonal.
        :param runs:    List containing locations of connected diagonal runs
        :param query:   Query sequence
        :param db:      Database sequence
        :param w:       Offset with resepect to diagonal
        :return:        PairAlignment object containing the local alignment
        """
        matrix = [[Pair(0) for j in range(len(db) + 1)]
                  for i in range(len(query) + 1)]

        highest = 0, 0, 0

        for i in range(1, len(matrix)):
            for j in range(1, len(matrix[0])):
                bounded = False
                for r in runs:
                    index = i-j
                    if r-w <= index <= r+w:
                        bounded = True

                    if bounded:
                        scores = [(matrix[i][j - 1] - self.gap_open), (matrix[i - 1][j] - self.gap_open), 0]

                        tentative_match = matrix[i - 1][j - 1].val
                        tentative_match += self.score_matrix.lookup(query[i - 1], db[j - 1])

                        scores.insert(0, tentative_match)

                        matrix[i][j].val = max(scores)
                        matrix[i][j].pred = scores.index(matrix[i][j].val)
                        if matrix[i][j].val > highest[0]:
                            highest = matrix[i][j].val, i, j

        i = highest[1]
        j = highest[2]

        aligned_s1 = ""
        aligned_s2 = ""

        score = 0

        while matrix[i][j].val > 0:
            if matrix[i][j].pred == 0:
                i -= 1
                j -= 1

                aligned_s1 += query[i]
                aligned_s2 += db[j]

                score += self.score_matrix.lookup(query[i], db[j])

            elif matrix[i][j].pred == 1:
                j -= 1

                aligned_s1 += "-"
                aligned_s2 += db[j]

                score -= self.gap_open

            else:
                i -= 1

                aligned_s1 += query[i]
                aligned_s2 += "-"

                score -= self.gap_open

        return PairAlignment(aligned_s1[::-1], aligned_s2[::-1], score,
                             s1_name=query.name, s2_name=db.name)


    def fasta(self, query_seq, db_seq, ktup, w):
        """
        FASTA alignment algorithm.
        :param query_seq:   Query sequence
        :param db_seq:      List of database sequences.
        :param ktup:        Length of k-words to search for.
        :param w:           Offset for banded Smith-Waterman algorithm.
        :return:            List of PairAlignment objects
                                containing local alignments.
        """
        chain_n = []
        init_n = []
        alignments = []
        query_kwords = {}

        # Generates neighbourhood of k-words
        for i in range(len(query_seq) - (ktup-1)):
            if query_seq[i:i+ktup] not in query_kwords:
                query_kwords[query_seq[i:i+ktup]] = []
            query_kwords[query_seq[i:i+ktup]].append(i)

        for other in db_seq:
            G = Graph()
            other_kwords = {}
            common_kwords = {}
            diag_starts = {}

            # Generates k-words for database sequence
            for i in range(len(other) - (ktup-1)):
                if other[i:i + ktup] not in other_kwords:
                    other_kwords[other[i:i + ktup]] = []
                other_kwords[other[i:i + ktup]].append(i)

            for word, rows in query_kwords.items():
                if word in other_kwords:
                    columns = other_kwords[word]
                    for i in rows:
                        for j in columns:
                            if i-j not in common_kwords:
                                common_kwords[i-j] = 0
                                diag_starts[i-j] = (i, j)
                            common_kwords[i-j] += self.match
            ordered_diags = sorted(common_kwords,
                                   key=common_kwords.get,
                                   reverse=True)[0:len(common_kwords) if
                                   len(common_kwords) <= 10 else 10]
            ordered_diags = {d: ordered_diags for i, d in enumerate(ordered_diags)}

            for i, diag in enumerate(ordered_diags):
                i = diag_starts[diag][0]
                j = diag_starts[diag][1]
                i_end = i + common_kwords[diag]
                j_end = j + common_kwords[diag]
                common_kwords[diag] = 0
                while i < min(i_end, len(query_seq)) and j < min(j_end, len(other[j])):
                    common_kwords[diag] += self.score_matrix.lookup(query_seq[i], other[j])
                    i += 1
                    j += 1
                ordered_diags[diag] = common_kwords[diag]
                G.add_vertex(diag, ordered_diags[diag])

            for diag, score in ordered_diags.items():
                for other_diag in ordered_diags:
                    if diag != other_diag:
                        i = diag_starts[other_diag][0] - (diag_starts[diag][0] + score)
                        j = diag_starts[other_diag][1] - (diag_starts[diag][1] + score)
                        if i >= 0 and j >= 0:
                            G.add_edge(diag, other_diag, -((i+j)*self.gap_open)//2)

            chain = G.chaining()
            chain_n.append(chain)
            init_n.append([ordered_diags[i] for i in chain])

        # PUT CODE HERE TO THRESHOLD INIT_N

        for i, run_list in enumerate(chain_n):

            alignments.append(self._banded_sw(
                                              run_list, query_seq, db_seq[i], w)
                              )
        return alignments
