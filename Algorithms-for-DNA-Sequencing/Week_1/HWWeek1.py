from typing import Tuple, List

class HWWeek1():
    """Class to solve the homework Week 1 of the course
    """
    def __init__(self, filename):
        self.filename = filename
    
    def read_fastq(self) -> Tuple[List, List]:
        """Reads a fastq file and returns the sequences and the qualities

        Returns
        -------
        Tuple[List, List]
            Tuple with the sequences and qualities in the fastq file
        """
        sequences = []
        qualities = []
        with open(self.filename) as f:
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
    
    def read_genome(self) -> str:
        """Reads a genome in a FASTA file

        Returns
        -------
        str
            genome
        """
        genome = ''
        with open(self.filename, 'r') as f:
            for line in f:
                # ignore header line with genome information
                if line[0] != '>':
                    genome += line.rstrip()
        return genome
    
    def reverse_complement(self, string: str):
        """Returns the reverse complement of a DNA string

        Parameters
        ----------
        s : str
            string of DNA
        """
        complements = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "N": "N"
        }
        t = ""
        for base in string:
            t = complements[base] + t # preppend to reverse string
        
        return t
    
    def phred_33_to_q(self, char: str) -> int:
        """_summary_

        Parameters
        ----------
        char : str
            _description_

        Returns
        -------
        int
            _description_
        """
        return ord(char) - 33

    def create_hist(self, qualities: list) -> list:
        """_summary_

        Parameters
        ----------
        qualities : list
            _description_

        Returns
        -------
        list
            _description_
        """
        hist = [0] * 50
        for qual in qualities:
            for phred in qual:
                q = HWWeek1.phred_33_to_q(phred)
                hist[q] += 1
        
        return hist
    
    def find_GC_by_pos(self, reads: list) -> list:
        """_summary_

        Parameters
        ----------
        reads : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        gc = [0] * 100
        totals = [0] * 100
        
        for read in reads:
            for i in range(len(read)):
                if read[i] == "C" or read[i] == "G":
                    gc[i] += 1
                totals[i] += 1
        
        for i in range(len(gc)):
            if totals[i] != 0:
                gc[i] /= float(totals[i])
                
        return gc

    def naive_with_rc(self, pattern: str, text: str) -> int:
        """_summary_

        Parameters
        ----------
        pattern : str
            _description_
        text : str
            _description_

        Returns
        -------
        int
            _description_
        """
        occurrences = []
        for i in range(len(text) - len(pattern) + 1):  # loop over alignments
            for seq in (pattern, HWWeek1.reverse_complement(self, pattern)):
                match = True
                for j in range(len(seq)):  # loop over characters
                    if text[i+j] != seq[j] :  # compare characters
                        match = False
                        break
                if match:
                    occurrences.append(i)
                    break    
    
        return occurrences

    def naive_2mm(self, pattern: str, text: str) -> list:
        """_summary_

        Parameters
        ----------
        pattern : str
            _description_
        text : str
            _description_

        Returns
        -------
        list
            _description_
        """
        occurrences = []
        for i in range(len(text) - len(pattern) + 1):
            mismatches = 0 # loop over alignments
            for j in range(len(pattern)):  # loop over characters
                if text[i+j] != pattern[j]:  # compare characters
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                occurrences.append(i)  # all chars matched; record
        return occurrences