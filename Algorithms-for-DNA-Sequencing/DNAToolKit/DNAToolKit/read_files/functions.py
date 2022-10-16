from typing import Tuple, List

def read_fastq(filename: str) -> Tuple[List, List]:
    """Reads a fastq file and returns the sequences and the qualities

    Parameters
    ----------
    filename: str
        filepath to a fastq file

    Returns
    -------
    Tuple[List, List]
        Tuple with the sequences and qualities in the fastq file
    """
    sequences = []
    qualities = []
    with open(filename) as f:
        while True:
            f.readline()  # skip name line
            seq = f.readline().rstrip()  # read base sequence
            f.readline()  # skip placeholder line
            qual = f.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def read_fasta(filename: str) -> str:
    """Reads a fasta file and stores the genome in a string

    Parameters
    ----------
    filename : str
        File path to the fasta file

    Returns
    -------
    str
        genome
    """
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
