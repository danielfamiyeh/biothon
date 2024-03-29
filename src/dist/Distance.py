class Distance:
    """
    Class with static methods for distance metrics on pairs of strings.
    """
    @staticmethod
    def hamming(str1, str2):
        """
        Hamming distance metric.
        :param str1:
        :param str2:
        :return:
        """
        dist = 0

        # Distances must be the same length for hamming distance
        if len(str1) != len(str2):
            print("String lengths must be equal for hamming distance.\nPlease use edit distance.")
        else:
            # Iterate over length of both strings
            for i in range(len(str1)):
                # Increment counter if indel is found
                if str1[i] != str2[i]:
                    dist += 1
        return dist

    @staticmethod
    def edit(str1, str2):
        """
        Levenshtein (edit) distance metric.
        :param str1:
        :param str2:
        :return:
        """
        rows = len(str1) + 1
        cols = len(str2) + 1
        matrix = []

        # Construct distance matrix
        for i in range(rows):
            matrix.append([])
            for j in range(cols):
                matrix[i].append(0)
                if j == 0:
                    matrix[i][j] = i
                elif i == 0:
                    matrix[i][j] = j

        # Traverse matrix
        for i in range(1, rows):
            for j in range(1, cols):
                # Get minimum number of operations per character
                matrix[i][j] = min(matrix[i-1][j-1] + (1 if str1[i-1] != str2[j-1] else 0),
                                   matrix[i-1][j] + 1,
                                   matrix[i][j-1] + 1)

        # Return count at last letter of both strings
        return matrix[rows-1][cols-1]

