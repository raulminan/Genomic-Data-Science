import regex as re
from typing import Dict, Tuple


class DNAToolKit():
    """This class contains the functions necessary to answer the questions in the 
    final exam for the course
    """
   
    def __init__(self, filepath: str) -> Dict[str, str]:
        """Opens a FASTA file and saves it in a dictionary, where the keys are
        the sequence identifiers and the values are the DNA sequences

        Parameters
        ----------
        filepath : str
            filepath to the FASTA file

        Returns
        -------
        dict
            Dictionary with sequence identifiers and the corresponding sequences
        """
        self.filepath = filepath
        self.dict = {}

        with open(self.filepath) as f:
            for line in f:
                line = line.rstrip()
                if line[0] == ">":
                    identifier = line
                    self.dict[identifier] = ""
                else:
                    self.dict[identifier] += line
    
    def record_count(self) -> int:
        """Counts the number of records in a FASTA file

        Returns
        -------
        int
            number of records
        """
        n = len(self.dict)
        return n
    
    def seq_lenght(self) -> Tuple[Dict, Dict]:
        """Obtains the sequences with the max and min lengths, their length and their
        identifier

        Returns
        -------
        Tuple[Dict, Dict]
            Tuple with the data for the sortest and longest sequence:

            Tuple['Shortest Sequence':{"identifier": identifier,
                                     "sequence": sequence,
                                     "length": length},

                 'Longest Sequence':{"identifier": identifier,
                                    "sequence": sequence,
                                    "length": length}]
        """
        max_len = 0
        min_len = 1e100

        max_seq = ""
        min_seq = ""

        longest_seq = {}
        shortest_seq = {}

        min_len_id = 0
        max_len_id = 0

        for id, seq in self.dict.items():
            if len(seq) > max_len:
                max_len = len(seq)
                max_len_id = id
                max_seq = seq
            if len(seq) < len(seq):
                min_len = len(seq)
                min_len_id = id
                min_seq = seq
        
        longest_seq["identifier"] = max_len_id
        longest_seq["sequence"] = max_seq
        longest_seq["length"] = max_len

        shortest_seq["identifier"]  = min_len_id
        shortest_seq["sequence"] = min_seq
        shortest_seq["length"] = min_len
    
        return shortest_seq, longest_seq
    
    def orf_finder(self, sequence: str, frame: int) -> Tuple[str, int, int]:
        """"Finds all open reading frames (ORFs) present in each sequence of the
        FASTA file

        Parameters
        ----------
        sequence : str
            a dna sequence to find the ORFs in
        frame : int
            reading grame to fing ORFs. Must be 1, 2 or 3

        Returns
        -------
        Tuple[str, int, int]
            Tuple with the longest ORF, its legnth and its position

        Raises
        ------
        ValueError
            raises ValueError if frame isn't 1, 2 or 3
        """
        # TODO fix bug
        if frame not in {1, 2, 3}:
            raise ValueError("Frame must be 1, 2 or 3")

        start = "ATG"
        stop = ["TAA", "TAG", "TGA"]

        max_orf_len = 0
        max_orf_start = 0
        max_orf = ""
        temp_max_orf = "" # stores temporary orfs
        sequence = sequence[frame - 1:]

        for i in range(0, len(sequence)-3, 3):
            if sequence[i:i+3] == start:
                for j in range(i, len(sequence)-3, 3):
                    codon = sequence[j:j+3]
                    temp_max_orf += codon
                    orf_len = len(sequence[i:j+3])
                    if codon in stop:
                        if orf_len > max_orf_len:
                            max_orf_len = orf_len
                            max_orf_start = i + frame # add the bases removed at sequence = sequence[frame-1]
                            max_orf = temp_max_orf
                            temp_max_orf = ""
                        break
        
        if max_orf == "":
            pass
        else:
            return max_orf, max_orf_len, max_orf_start
        

    def repeats(self, sequence: str, n: int) -> dict:
        """Given a length n, identifies all repeats of length n in all sequences of a 
        FASTA file


        Parameters
        ----------
        sequence : str
            sequence of dna to find repeats in 
        n : int
            length of the repeats to find

        Returns
        -------
        dict
            dictionary with the repeats and the number of times it appears
        """
        repeats = {}

        for i in range(0, len(sequence)-n+1):
            repeat = sequence[i:i+n]
            matches = re.findall(repeat, sequence, overlapped=True)
            repeats[repeat] = len(matches)

        return repeats
