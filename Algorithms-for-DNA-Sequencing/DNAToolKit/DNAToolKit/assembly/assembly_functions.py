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
    """Creates an overlap map from a list of DNA strings

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
