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
        """Turn Phred+33 ASCII-encoded quality into Q

        Parameters
        ----------
        char : str
            ASCII encoded quality

        Returns
        -------
        int
            Q
        """
        return ord(char) - 33

    def q_to_phred_33(self, q: int) -> str:
        """Turns Q value into a Phred+33 ASCII-encoded quality char

        Parameters
        ----------
        q : int
            Q value

        Returns
        -------
        str
            Phred+33 ASCII-encoded representation of the Q value
        """
        return chr(q+33)

    def create_hist(self, qualities: list) -> list:
        """Given a list of qualities, builds a histogram with
        the frequency of each quality in qualtities

        Parameters
        ----------
        qualities : list
            list of qualities

        Returns
        -------
        list
            list of the frequency of each quality
        """
        hist = [0] * 50
        for qual in qualities:
            for phred in qual:
                q = HWWeek1.phred_33_to_q(phred)
                hist[q] += 1
        
        return hist
    
    def find_GC_by_pos(self, reads: list) -> list:
        """Given a list of reads, returns the G/C content
        at each position

        Parameters
        ----------
        reads : list
            list of reads

        Returns
        -------
        list
            list with the G/C content at each position
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

    def naive_with_rc(self, pattern: str, text: str) -> list:
        """Naive matching betweem a pattern and a text in both strands
        of a DNA string

        Parameters
        ----------
        pattern : str
            dna string to find in text
        text : str
            dna string where pattern wil be found

        Returns
        -------
        list
            list of occurences of the matches
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
        """Naive matching with up to 2 mismatches

        Parameters
        ----------
        pattern : str
            dna string to find in text
        text : str
            dna string where pattern wil be found

        Returns
        -------
        list
            list of occurences of the matches
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
        