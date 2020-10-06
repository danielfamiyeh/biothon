def get_composition(string, k, **kwargs):
    """
    Returns k-mer composition of a given string
    :param kwargs:
        sorted: (Boolean): Sorts composition lexicographically if True,
                                else does not
    :param string: String to get composition of
    :param k: length of k-mers in composition
    :return: k-mer composition of string
    """
    sorted = kwargs.get("sorted", False)
    kmer_comp = [string[i: i + k] for i in range(len(string) - (k - 1))]
    if sorted:
        kmer_comp.sort()
    return kmer_comp
