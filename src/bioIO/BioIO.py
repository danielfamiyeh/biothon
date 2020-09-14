from src.seq.Seq import *
from src.align.msa import MultiAligner
from src.align.msa import MSAProfile
from src.align.pair import PairAlignment
from copy import deepcopy

def _write_fasta(seq, file):
    """
    Writes sequences or a pairwise alignment to a .fasta file.
    :param seq:     Sequence/list of sequences/pairwise alignment to write to file.
    :param file:    Filename
    :return:        None
    """
    if isinstance(seq, Seq):
        file.write(f">{str(seq).replace(':', '')}\n")
    elif isinstance(seq, list):
        for s in seq:
            file.write(f">{str(s).replace(':', '')}\n")
            if len(s) != 60:
                file.write("\n")


def _write_clustal(seq, file):
    """
    Writes a multiple sequence alignment to a .clustal file
    :param seq:     MSA to write to file in the form of a MultiAligner or
                    MSAProfile object
    :param file:    Filename
    :return:        None
    """
    if isinstance(seq, MultiAligner):
        file.write(f"{seq.profile}")
        return
    file.write(f"{seq}")


# Validate seq type
class BioIO:
    """
    Biothon I/O Class.
        Contains static methods for file I/O.
    """
    @staticmethod
    def read(fn, seq_type):
        """
        Reads a .fasta or .clustal file into a sequence or MSPofile respectively
        :param fn:          Filename
        :param seq_type:    SeqType denoted by data in file.
        :return:            None
        """
        # Check that param: fn is str and of correct format
        if isinstance(fn, str):
            if "." in fn:
                name, ext = fn.split(".")
                # If file is a fasta file
                if ext == "fasta":
                    # Initialise empty sequence list
                    seqs = []
                    with open(fn, "r") as fasta_file:
                        seq_name, seq_desc, seq = "", "", ""
                        # Enumerate lines in file
                        for i, line in enumerate(fasta_file):
                            # Line is the header line
                            if line.startswith(">"):
                                if i > 0:
                                    seqs.append(Seq(seq, seq_type, name=seq_name, desc=seq_desc))
                                    seq_name, seq_desc, seq = "", "", ""
                                if "|" in line:
                                    pipe_index = line.find("|")
                                    seq_name = line[1:pipe_index-(1 if line[pipe_index-1] == " " else 0)]
                                    seq_desc = line[pipe_index + (2 if line[pipe_index + 1] == " " else 1):]
                                else:
                                    seq_name = line
                            else:
                                seq += line.strip()
                        if len(seqs) == 0 and seq != "":
                            seqs = Seq(seq, seq_type, name=seq_name, desc=seq_desc)
                    return seqs

                elif ext == "clustal":
                    with open(fn, "r") as clustal_file:
                        j = 0
                        names = []
                        seqs = []
                        for i, line in enumerate(clustal_file):
                            if i > 1 and line[0] != '\n':
                                after_name = line.find("\t")
                                name = line[0: after_name].rstrip()
                                partial_seq = line[after_name:].strip()
                                partial_seq = partial_seq[0: partial_seq.find(" ")]

                                if len(names) < j + 1:
                                    names.append(name)
                                    seqs.append(partial_seq)
                                else:
                                    seqs[j] += partial_seq
                                j += 1
                            else:
                                j = 0

                        profile = MSAProfile(Seq("", seq_type, name=""), AminoScoreMatrix(AminoScoreType.BLOSUM62))
                        profile.seq_list = deepcopy([Seq(seqs[i], seq_type, name=names[i]) for i in range(len(seqs))])
                        profile.reprofile()
                        return profile
                else:
                    raise ValueError("File extension must be fasta or clustal.")
            else:
                raise ValueError("Param fn must have a file extension.")
        else:
            raise TypeError("Param fn must be of type string.")

    @staticmethod
    def write(seq, fn):
        """
        Static method to write a Seq/list of Seq pbkects/pariwise or MSA to file
        :param seq: Object to write to file
        :param fn:  Filename
        :return:    None
        """
        ext = ""
        if isinstance(fn, str):
            if isinstance(seq, Seq) or isinstance(seq, PairAlignment):
                ext += ".fasta"
            elif isinstance(seq, list):
                for i, s in enumerate(seq):
                    if not isinstance(s, Seq):
                        raise TypeError("Item {i} in list is type {type(s)}.\n"
                                        "All items must be of type Seq")
                ext += ".fasta"
            elif isinstance(seq, MultiAligner) or isinstance(seq, MSAProfile):
                ext += ".clustal"
            else:
                raise TypeError("Param seq must be of type Seq,"
                                "a list of type Seq or an MSAProfile/MultiAligner.")
            with open(f"{fn if '.' not in fn else fn[0:fn.index('.')]}{ext}", "w") as file:
                if ext == ".fasta":
                    _write_fasta(seq, file)
                else:
                    _write_clustal(seq, file)

        else:
            raise TypeError("Param fn must be a string.")



