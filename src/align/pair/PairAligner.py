from math import ceil
from src.align.pair import PairAlignment
from src.seq.Seq import *
from collections import namedtuple
from src.graph import Graph

class Pair:
    def __init__(self, val):
        """
        Pair object constructor, used to represent entry in matrix
        :param val: Value of entry in matrix.
        """
        self.val = val
        self.pred = None    # Predecessor object reference

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
        :return:                Returns PairAlignment object
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
                    # Add indels
                    scores = [(matrix[i][j-1]-self.gap_open),
                              (matrix[i-1][j]-self.gap_open)]

                    # Add tentative match
                    tentative_match = matrix[i - 1][j - 1].val
                    tentative_match += self.match if s1[i - 1] == s2[j - 1] else -self.mismatch

                    # A match case is preferred over all other cases
                    scores.insert(0, tentative_match)

                    # Get argmax and set predecessor
                    matrix[i][j].val = max(scores)
                    matrix[i][j].pred = scores.index(matrix[i][j].val)

            # Start from bottom-right of matrix
            i = rows - 1
            j = cols - 1

            # Aligned strings
            aligned_s1 = ""
            aligned_s2 = ""

            # For affine-gap penalty
            k = 0

            # Iterate until start
            while i > 0 or j > 0:
                # If no more letters in query string
                if i == 0:
                    j -= 1
                    for count in range(j, -1, -1):
                        k += 1
                        j -= 1
                        # Add gaps to query string
                        aligned_s1 += "-"
                        # Add the rest of letters to db string
                        aligned_s2 += s2[count]
                        # Subtract gap penalty
                        score -= self.gap_open

                # If no more letters in db string
                elif j == 0:
                    i -= 1
                    for count in range(i, -1, -1):
                        k += 1
                        i -= 1

                        # Add rest of letters to query string
                        aligned_s1 += s1[count]
                        # Add gaps to db string
                        aligned_s2 += "-"
                        # Subtract gap penalty
                        score -= self.gap_open

                # If predecessor was a match
                elif matrix[i][j].pred == 0:
                    # Backtrack to upper-left corner neighbour
                    i -= 1
                    j -= 1
                    k = 0

                    # Append character to both sequences
                    aligned_s1 += s1[i]
                    aligned_s2 += s2[j]

                    # Increment score by based on value from score matrix
                    score += self.score_matrix.lookup(s1[i], s2[j])

                # If predecessor was a query indel
                elif matrix[i][j].pred == 1:
                    # Backtrack to left neighbour in alignment matrix
                    j -= 1
                    k += 1

                    # Add gap to query string
                    aligned_s1 += "-"
                    # Add character to db string
                    aligned_s2 += s2[j]

                    # Subtract gap penalty
                    score -= self.gap_open

                # If predecessor was a db indel
                else:
                    # Move to neighbour above in alignment matrix
                    i -= 1
                    k += 1

                    # Add character to query string
                    aligned_s1 += s1[i]

                    # Add gap to db string
                    aligned_s2 += "-"

                    # Subtract gap penalty
                    score -= self.gap_open

            # Return new pair alignment object with correct directionality
            return PairAlignment(aligned_s1[::-1], aligned_s2[::-1], score,
                                 s1_name=s1.name, s2_name=s2.name)
        else:
            # Throw type error
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

            # Alignment matrix rows and columns
            rows = len(s1) + 1
            cols = len(s2) + 1

            matrix = []

            # Generate alignment matrix with all zeroes
            for i in range(rows):
                matrix.append([Pair(0) for _ in range(cols)])

            # Iterate over alignment matrix rows
            for i in range(1, rows):
                # Iterate over alignment matrix columns
                for j in range(1, cols):
                    # And indel scores to score list
                    scores = [(matrix[i][j-1] - self.gap_open), (matrix[i-1][j] - self.gap_open), 0]

                    # Set tentative match to previous matches score
                    tentative_match = matrix[i-1][j-1].val
                    # Add score based on score matrix
                    tentative_match += self.match if s1[i-1] == s2[j-1] else -self.mismatch
                    # Add tentative match to first index of score list
                    scores.insert(0, tentative_match)

                    # Set value of value of currrent element in alignment matrix
                    # to argmax(scores)
                    matrix[i][j].val = max(scores)
                    # Set predecessor value to first index of argmax(scores)
                    matrix[i][j].pred = scores.index(matrix[i][j].val)

            # Highest value in alignment matrix
            highest = 0
            # Row of max value in alignment matrix
            max_row = 0
            # Column of max value in alignment matrix
            max_col = 0

            # Local alignment strings
            s1_local = ""
            s2_local = ""

            # TODO: For affine gap penalty?
            k = 0

            # Iterate over alignment matrix rows
            for i, row in enumerate(matrix):
                # Iterate over alignment matrix columns
                for j, col in enumerate(row):
                    # If highest is less than
                    # the current element
                    if highest < col.val:
                        # Set highest to current element's value
                        highest = col.val
                        # Set row and column
                        max_row = i
                        max_col = j

            # Iterate through matrix while current value > 0
            while matrix[max_row][max_col].val > 0:
                # If predecessor was a match
                if matrix[max_row][max_col].pred == 0:
                    # Move to upper left neighbour in alignment matrix
                    max_row -= 1
                    max_col -= 1
                    k = 0

                    # Add character to both local alignment strings
                    s1_local += s1[max_row]
                    s2_local += s2[max_col]

                    # Increment score by score matrix value
                    score += self.score_matrix.lookup(s1[max_row], s2[max_col])

                # If predecessor was a query indel
                elif matrix[max_row][max_col].pred == 1:
                    # Move to left neighbour in alignment matrix
                    max_col -= 1
                    k += 1

                    # Add gap to query string
                    s1_local += "-"
                    # Add character to db string
                    s2_local += s2[max_col]

                    # Substract gap penalty
                    score -= self.gap_open

                # Else predecessor was a db indel
                else:
                    # Move to neighbour above in alignment matrix
                    max_row -= 1
                    k += 1

                    # Add character to query string
                    s1_local += s1[max_row]
                    # Add gap to db string
                    s2_local += "-"

                    # Subtract gap penalty
                    score -= self.gap_open

            # Return PairAlignment with correct directionality
            return PairAlignment(s1_local[::-1], s2_local[::-1], score,
                                 s1_name=s1.name, s2_name=s2.name)
        else:
            raise TypeError("Params s1 and s2 must be of type Seq.")

    def __needle_score(self, s1, s2):
        """

        :param s1:              First query sequence
        :param s2:              Second query sequence
        :return:                Returns last line of Needleman-Wunsch alignment.
        """
        rows = len(s1) + 1
        cols = len(s2) + 1

        # Generate alignment matrix of height 2
        matrix = [[0 if i == 1 else j * -self.gap_open for j in range(cols)] for i in range(2)]

        for i in range(1, rows):
            matrix[1][0] = i * -self.gap_open
            for j in range(1, cols):

                # Add indels to score list
                scores = [(matrix[1][j - 1] - self.gap_open), (matrix[0][j] - self.gap_open)]
                # Set tentative match to previous matches' score
                tentative_match = matrix[0][j-1]
                # Increment value of tentative match by match value
                # TODO: Change to score matrix
                tentative_match += self.match if s1[i - 1] == s2[j - 1] else -self.mismatch

                # Insert tentative match into first index of scor elist
                scores.insert(0, tentative_match)
                # Set value of current entry to argmax(scores)
                matrix[1][j] = max(scores)

            # Move down alignment matrix
            matrix[0] = matrix[1].copy()
            matrix[1] = [0 for _ in range(cols)]

        # Return new first row
        return matrix[0]

    def _hirschberg_global(self, s1, s2, score):
        """
        Internal hirschberg algorithm that gets called recursively.
        :param s1: Sequence to align
        :param s2: Sequence to align
        :param score: Running total of score
        :return: 3-tuple containing two aligned strings and a score
        """

        # Aligned strings
        str1 = ""
        str2 = ""

        # If we are at beginning of query string
        if len(s1) == 0:
            # Iterate over characters left in db string
            for char in s2:
                # Add gaps to query string
                str1 += "-"
                # Add characters to db string
                str2 += char
                # Subtract gap penalty
                score -= self.gap_open

        # If we are at the beginning of db string
        elif len(s2) == 0:
            # Iterate over characters left in query string
            for char in s1:
                # Add characters to query string
                str1 += char
                # Add gaps to db string
                str2 += "-"
                # Subtract gap penalty
                score -= self.gap_open

        # If either string has one character left
        elif len(s1) == 1 or len(s2) == 1:
            # Get NW alignment of both
            t = self.needle(s1, s2)
            # Append to query string
            str1 += t.seq1
            # Append to db string
            str2 += t.seq2
            # Increment score by alignment score
            score += t.score

        # Else perform Hirschberg's algorithm
        else:
            # Get midpoint
            y_mid = ceil(len(s1) // 2)

            # Divide matrix into 2 halves about the midpoint
            top = self.__needle_score(s1[0:y_mid], s2)
            bot = self.__needle_score((s1[y_mid:])[::-1], s2[::-1])[::-1]
            # Sum element-wise scores
            sum_vec = [top[i] + bot[i] for i in range(len(s2))]
            # Calculate x midpoint
            x_mid = sum_vec.index(max(sum_vec))

            # Recursive splits
            t1 = self._hirschberg_global(Seq(s1[0:y_mid], s1.seq_type), Seq(s2[0: x_mid], s2.seq_type), score)
            t2 = self._hirschberg_global(Seq(s1[y_mid:], s1.seq_type), Seq(s2[x_mid:], s2.seq_type), score)

            # Append to query string
            str1 += t1[0] + t2[0]
            # Append to db string
            str2 += t1[1] + t2[1]
            # Increment score
            score += t1[2] + t2[2]

        return str1, str2, score

    def hirschberg_global(self, s1, s2):
        """
        Hirschberg's algorithm for global alignment.
        Computes an optimised Needleman-Wunsch alignment on a pair of strings.
        :param s1: Query string
        :param s2: DB string
        :return: PairAlignment object or throws TypeError
        """
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

        # Query string
        query = kwargs["query"]
        # Database string list
        db = kwargs["db"]
        # Length of k-words
        ktup = kwargs["ktup"]
        # Thresold value for k-word consideration
        threshold = kwargs["threshold"]
        # Cuttof point for high-scoring value
        cutoff = kwargs["cutoff"]
        # Minimum value for HSP consideration
        min_score = kwargs["min_score"]

        # Initialise highest-score to zero
        highest_score = 0

        # Initialise empty neighbourhood
        neighbourhood = {}

        # Initialise empty set of seed hits
        hits = {i: set() for i in range(len(db))}
        hsp = []

        # Initialise alpha based on the query strings' alphabet
        alpha = alphabet_dna if query.seq_type is SeqType.DNA else alphabet_rna \
            if query.seq_type is SeqType.RNA else alphabet_protein_1

        # Populating k-word neighbourhood
        # Iterate over query string wrt ktup
        for i in range(len(query) - ktup):
            # For every character in query's alphabet
            for c in alpha:
                # Initialise new k-word
                kword = c + query[i+1:ktup]
                # Initialise query substring
                query_sub = query[i:ktup]

                # Initialise number of matches to zero
                match = 0
                # Iterate over query substring
                for j in range(len(query_sub)):
                    # Increment match value based on score matrix
                    match += self.score_matrix.lookup(query_sub[j], kword[j])

                # If match score breaks threshold
                if match > threshold:
                    # If index is not a key in neighbourhood map
                    if i not in neighbourhood:
                        # Add an empty list if key i
                        neighbourhood[i] = []
                    # Add k-word to list in heighbourhood
                    neighbourhood[i].append(kword)

        # Getting seed hits
        # Enumerate database sequences
        for i, seq in enumerate(db):
            # Iterate over neighbourhood map
            for key, value in neighbourhood.items():
                # Enumerate every list in neighbourhood
                for j, neighbour in enumerate(value):
                    # Initialise character count
                    k = 0
                    while True:
                        # If character count greater than length of db string
                        if k > len(seq):
                            # Exit loop
                            break
                        # Initialise match index
                        match = seq[k: k+ktup].find(neighbour)
                        # If match between two sequences
                        if match != -1:
                            # Add indices of character match to seed hit
                            hits[i].add((key, j, k + match))
                        # Increment character count
                        k += 1

        # Initialise high scoring pair as a namedtuple.
        HSPTuple = namedtuple("HSP", [f"seq_label", "subseq", "score", "query_index", "seq_index"])

        # Getting high scoring pairs
        # Iterate over seed hit list
        for i, hit_list in hits.items():
            # If list is not empty
            if len(hit_list) > 0:
                # Enumerate seed hit list
                for j, seed_hit in enumerate(hit_list):
                    # Left/right 'panes' of search window
                    left = 0
                    right = 0

                    # TODO: Infinite loops
                    while True:
                        # Assign query and db subsequences
                        query_sub = query[seed_hit[0] - left: seed_hit[0] + ktup + right]
                        seq_sub = db[i][seed_hit[2] - left: seed_hit[2] + ktup + right]

                        # Initialise match score
                        match = 0
                        # Iterate over characters in sequences
                        for k in range(min(len(seq_sub), len(query_sub))):
                            # Increment match score using score matrix
                            match += self.score_matrix.lookup(query_sub[k], seq_sub[k])

                        # If match value signifies a high-scoring pair
                        if match >= min_score and match > highest_score:
                            hsp.append(HSPTuple((db[i].label if db[i].label != "" else i),
                                                seq_sub, match, seed_hit[0] - left,
                                                seed_hit[2] - left))

                            # Update highest_score attained
                            highest_score = match

                        # If score has dropped below cutoff break
                        elif match < highest_score - cutoff:
                            break

                        # If extensions result in a negative score break
                        elif min(seed_hit[0] - left, seed_hit[2] - left) < 0 and\
                            (seed_hit[2] + ktup + right > len(db[i]) or
                             seed_hit[0] + ktup + right > len(query)):
                            break

                        # Extend window
                        left += 1
                        right += 1

        # Initialise ordered high scoring pair list
        hsp = sorted(hsp, key=lambda x: x.score, reverse=True)
        # Return hsp list
        return hsp

    def _banded_sw(self, runs, query, db, w):
        """
        Banded Smith-Waterman algorithm. Performs local alignment that is
        constrained to a certain distance from a diagonal.
        :param runs:    List containing locations of connected diagonal runs
        :param query:   Query sequence
        :param db:      Database sequence
        :param w:       Offset with respect to diagonal
        :return:        PairAlignment object containing the local alignment
        """

        # Generates an alignment matrix with all zeroes
        matrix = [[Pair(0) for j in range(len(db) + 1)]
                  for i in range(len(query) + 1)]

        # Initialise highest-score with indices tuple to all zeroes
        highest = 0, 0, 0

        # Iterate over rows of alignment matrix
        for i in range(1, len(matrix)):
            # Iterate over columns of alignment matrix
            for j in range(1, len(matrix[0])):
                # Assume we are not bounded
                bounded = False
                # Iterate over all runs
                for r in runs:
                    # Get diagonal's index
                    index = i-j
                    # If index is within distance of diagonal
                    if r-w <= index <= r+w:
                        # We are bounded
                        bounded = True

                    # If are bounded then perform Smith-Waterman alignment
                    if bounded:
                        # Add scores to score list, with zero so alignment can begin
                        # anywhere in the alignment matrix
                        scores = [(matrix[i][j - 1] - self.gap_open), (matrix[i - 1][j] - self.gap_open), 0]

                        # Initialise tentative_match to the previous match score value
                        tentative_match = matrix[i - 1][j - 1].val
                        # Increment according to value in score matrix
                        tentative_match += self.score_matrix.lookup(query[i - 1], db[j - 1])

                        # Insert new tentative match value to first index of
                        # score list
                        scores.insert(0, tentative_match)

                        # Set current value to argmax(scores)
                        matrix[i][j].val = max(scores)
                        # Set predecessor to index of first occurence of
                        # argmax(scores)
                        matrix[i][j].pred = scores.index(matrix[i][j].val)

                        # If current value is greater than the highest attained
                        if matrix[i][j].val > highest[0]:
                            # Update highest with new value and indices
                            highest = matrix[i][j].val, i, j

        # TODO: Why did I do this?
        i = highest[1]
        j = highest[2]

        # Initialise strings
        aligned_s1 = ""
        aligned_s2 = ""

        # Initialise score to zero
        score = 0

        # Backtrack over alignment matrix
        while matrix[i][j].val > 0:
            # If predecessor was a match
            if matrix[i][j].pred == 0:
                # Move to upper-left neighbour in alignment matrix
                i -= 1
                j -= 1

                # Add characters to both alignment strings
                aligned_s1 += query[i]
                aligned_s2 += db[j]

                # Increment score according to scor ematrix
                score += self.score_matrix.lookup(query[i], db[j])

            # If predecessor was a query indel
            elif matrix[i][j].pred == 1:
                # Move to left neighbour in alignment matrix
                j -= 1

                # Add gap to query aligned sequence
                aligned_s1 += "-"
                # Add character to db aligned sequence
                aligned_s2 += db[j]

                # Subtract score by gap penalty
                score -= self.gap_open

            # If predecessor was a db indel
            else:
                # Move to neighbour above in alignment matrix
                i -= 1

                # Add character to query string
                aligned_s1 += query[i]
                # Add gap to db string
                aligned_s2 += "-"

                # Subtract score by gap penalty
                score -= self.gap_open

        # Return pair alignment object with directionality corrected
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
        # Iterate over query sequence length wrt ktup
        for i in range(len(query_seq) - (ktup-1)):
            # If substring is not a key in query k-word map
            if query_seq[i:i+ktup] not in query_kwords:
                # Add a new empty list with substring as a key
                query_kwords[query_seq[i:i+ktup]] = []
            # Add index to list
            query_kwords[query_seq[i:i+ktup]].append(i)

        # Iterate over db sequences
        for other in db_seq:
            # Initialise empty graph
            G = Graph()
            # Initialise db k-word map
            other_kwords = {}
            # Initialise common k-word map
            common_kwords = {}
            # Initalise diagonal starts map
            diag_starts = {}

            # Generates k-words for database sequence
            # Iterate over db sequence length wrt ktp
            for i in range(len(other) - (ktup-1)):
                # If substring is not a key in other kwords map
                if other[i:i + ktup] not in other_kwords:
                    # Add substr as key with empty list as value
                    other_kwords[other[i:i + ktup]] = []
                # Add index to list
                other_kwords[other[i:i + ktup]].append(i)

            # Iterate over query k-words
            for word, rows in query_kwords.items():
                # If word is also present in db sequence
                if word in other_kwords:
                    # Initialise number of columns
                    # wrt db sequence k-word count
                    columns = other_kwords[word]
                    # Iterate over query rows
                    for i in rows:
                        # Iterate over db columns
                        for j in columns:
                            # If diag index is int present in
                            # common k-word map
                            if i-j not in common_kwords:
                                # Add it to map with value zero
                                common_kwords[i-j] = 0
                                # Add index as tuple to diag_starts map
                                diag_starts[i-j] = (i, j)
                            # Increment common k-word score
                            common_kwords[i-j] += self.match
            # Order diagonals by score and store in list
            ordered_diags = sorted(common_kwords,
                                   key=common_kwords.get,
                                   reverse=True)[0:len(common_kwords) if
                                   len(common_kwords) <= 10 else 10]

            # Dict comprehension by enumerating ordered_diags list
            ordered_diags = {d: ordered_diags for i, d in enumerate(ordered_diags)}

            # Enumerate ordered diagnoal dict
            for i, diag in enumerate(ordered_diags):
                # Set indices
                i = diag_starts[diag][0]
                j = diag_starts[diag][1]
                # Set ends
                i_end = i + common_kwords[diag]
                j_end = j + common_kwords[diag]

                # Initialise score counter
                common_kwords[diag] = 0

                # While within bounds of both sequences
                while i < min(i_end, len(query_seq)) and j < min(j_end, len(other[j])):
                    # Increment counter via score matrix alues
                    common_kwords[diag] += self.score_matrix.lookup(query_seq[i], other[j])
                    # Increment indices
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

            # Solve chaining problem and store result
            chain = G.chaining()
            # Add chain value to chain_n list
            chain_n.append(chain)
            # Add scores to init_n list
            init_n.append([ordered_diags[i] for i in chain])

        #TODO: PUT CODE HERE TO THRESHOLD INIT_N

        # Enumerate chain list
        for i, run_list in enumerate(chain_n):
            # Append alignments to alignment list
            alignments.append(self._banded_sw(
                                              run_list, query_seq, db_seq[i], w)
                              )
            # Return alignment list
        return alignments
