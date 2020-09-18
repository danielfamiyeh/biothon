import os

# Aboslute path of matrix module
dir = os.path.dirname(os.path.abspath(__file__))


def load(type_enum):
    """
    Load matrix from file.
    :param type_enum: Substitution matrix enum
    :return: Matrix dict
    """
    filename = str(type_enum).split(".")[1].lower()
    folder = ''.join([c for c in filename if not c.isdigit()])

    with open(f"{dir}\\..\\matrix\\{folder}\\{filename}.txt", "r") as matrix_file:
        keys = []
        mat = {}
        for i, line in enumerate(matrix_file):
            stripped = line.lstrip(' ')
            stripped = stripped[stripped.find(' '):]
            if i > 0:
                if i == 1:
                    keys = stripped.lstrip(' ').rstrip('\n').split('  ')
                    mat = {k: {k: 0 for k in keys} for k in keys}

                else:
                    single_spaced = stripped.replace("  ", " ").lstrip().rstrip('\n')
                    scores = single_spaced.split(' ')[1:]
                    for j in range(len(scores)):
                        mat[keys[i-2]][keys[j]] = scores[j]
        return mat
