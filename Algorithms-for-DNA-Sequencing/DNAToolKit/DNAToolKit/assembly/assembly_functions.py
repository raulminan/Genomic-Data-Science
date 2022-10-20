"""functions for genome assembly"""
from itertools import permutations
from typing import Tuple

def overlap(a: str, b: str, min_length: int = 3):
    """Returns the length of the longest overlap between the prefix of one
    string and the suffix of another string

    Parameters
    ----------
    a : str
        first string
    b : str
        second string
    min_length : int, optional
        minimum length of overlap, by default 3
    """
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        
        if b.startswith(a[start:]): # verify that the prefix of b is equal to the suffix of a
            return len(a) - start # length of overlap
        
        start += 1
    
def naive_overlap_map(reads: list[str], min_length: int) -> dict[Tuple[str, str], int]:
    """Creates an overlap map from a list of DNA strings followin a naive approach

    Parameters
    ----------
    reads : list[str]
        list of dna reads
    min_length : int
        minimum length of overlapped to be accounted for

    Returns
    -------
    dict[Tuple[str, str], int]
        dictionary with the map
        Example: {("AGG", "GGT"): 2} 
        AGG suffix and GGT prefix overlap by 2 characters 
    """
    overlaps = {}
    for a, b in permutations(reads, 2):
        overlap_length = overlap(a, b, min_length)
        if overlap_length > 0:
            overlaps[(a, b)] = overlap_length
    return overlaps

def overlap_map(reads: list[str], min_length: int) -> dict[Tuple[str, str], int]:
    """Creates an overlap map from a list of DNA strings. 

    Parameters
    ----------
    reads : list[str]
        list of dna reads
    min_length : int
        minimum length of overlapped to be accounted for

    Returns
    -------
    dict[Tuple[str, str], int]
        dictionary with the map
        Example: {("AGG", "GGT"): 2} 
        AGG suffix and GGT prefix overlap by 2 characters 
    """
    genome = "".join(reads)
    sets = {}
    overlaps = {}
    
    # let every kmer in the dataset have an associated set() object
    for i in range(len(genome), min_length + 1):
        kmer = genome[i:i+min_length]
        sets[kmer] = set()

    # for every kmer in a read, add the read to the set object that
    # corresponds to that kmer
    for read in reads:
        for i in range(len(read) - min_length + 1):
            kmer = read[i:i+min_length]
            try:
                sets[kmer].add(read)
            except KeyError:
                sets[kmer] = set()
                sets[kmer].add(read)
    
    # for each read "a", find all overlaps involving a suffix of "a"
    # take a's length-k suffix, find all reads containing that kmer,
    # then call overlap for each read
    for read in reads:
        suffix = read[-min_length:]
        kmers = list(sets[suffix])
        for kmer in kmers:
            if read == kmer:
                continue # don't overlap reads with themselves
            overlap_length = overlap(read, kmer, min_length)
            overlaps[(read, kmer)] = overlap_length

    return overlaps