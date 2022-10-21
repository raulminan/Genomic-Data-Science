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
    overlap_map = {}
    overlap_graph = {}
    
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
            if read != kmer:
                overlap_length = overlap(read, kmer, min_length)
                overlap_map[(read, kmer)] = overlap_length
                overlap_graph[read] = kmer

    return overlap_map, overlap_graph

def scs(ss: list) -> Tuple[str, list[str]]:
    """Computes the shortest common superstring (scs) of a set of strings
    The scs is the shortest string that contains all other strings as substrings 

    Parameters
    ----------
    ss : list
        set of strings

    Returns
    -------
    Tuple[str, list[str]]
        str
            first scs found
        list[str]
            list of all possible scs
    """
    scs = None
    scs_list = set()
    for permutation in permutations(ss):
        sup = permutation[0]
        for i in range(len(ss) - 1):
            overlap_length = overlap(permutation[i], permutation[i+1], 1)
            sup += permutation[i+1][overlap_length:] # append the part that doesn't overlap
    
        if scs is None or len(sup) < len(scs):
            scs = sup
            scs_list.add(scs)
        
        if len(sup) == len(scs):
            scs_list.add(sup)

    return scs, sorted(list(scs_list))

def pick_maximal_overlap(reads: list, min_overlap: int) -> Tuple[str, str, int]:
    """Returns the two strings with the maximum overlap and the overlap
    from a set of strings

    Parameters
    ----------
    reads : list
        set of dna reads
    min_overlap : int
        minium overlap to count

    Returns
    -------
    Tuple[str, str, int]
        two strings with the max overlap length and the legnth
    """
    read_a, read_b = None, None
    best_overlap_length = 0
    for a, b in permutations(reads, 2):
        overlap_length = overlap(a, b, min_length=min_overlap)
        if overlap_length > best_overlap_length:
            read_a, read_b = a, b
            best_overlap_length = overlap_length
    
    return read_a, read_b, best_overlap_length

def greedy_scs(reads: list, min_length: int) -> str:
    """Returns the shortest common superstring of a set of reads using a
    greedy approach

    Parameters
    ----------
    reads : list
        list of reads
    min_length : int
        min length of overlap allowed

    Returns
    -------
    str
        shortest common supersring of the set of reads
    """
    read_a, read_b, overlap_length = pick_maximal_overlap(reads, min_length)
    while overlap_length > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[overlap_length:]) # append suffix of b to a
        read_a, read_b, overlap_length = pick_maximal_overlap(reads, min_length)
    
    return "".join(reads) # join remaiming reads


def de_bruijn(read: str, k: int) -> dict:
    """Constructs a De Buijn graph from a string

    Parameters
    ----------
    read : str
        string that will be used for creating the De Bruijn graph
    k : int
        size of k-mer to construct the graph. The graph will contain one edge
        per kmer and one node for every distinct (k-1)-mer.

    Returns
    -------
    dict
        _description_
    """
    edges = []
    nodes = set()

    for i in range(len(read) - k + 1):
        edges.append((read[i:i+k-1], read[i+1:i+k]))
        nodes.add(read[i:i+k-1])
        nodes.add(read[i+1:i+k])
    return nodes, edges
