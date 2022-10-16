#!/usr/bin/env python
import os

from HomeWorkWeek2 import HomeWorkWeek2
from bm_preproc import BoyerMoore


if __name__ == "__main__":
    hw = HomeWorkWeek2
    FILENAME = os.path.join(os.path.dirname(__file__), "data", "chr1.GRCh38.excerpt.fasta")
    TEXT = hw.read_genome(FILENAME)

    # Q1, Q2
    # How many alignments and character comparisons does the naive exact matching 
    # algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
    # (derived from human Alu sequences) to the excerpt of human chromosome 1?  
    # (Don't consider reverse complements.)
    
    PATTERN_1 = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
    occurrences, num_aligments, num_character_comparisons =\
        hw.naive_with_counts(PATTERN_1, TEXT)
        
    print(f"The naive exact matching algorithm tries {num_aligments} alignments and "
        f"{num_character_comparisons} character comparisons")
    
    # Q3
    # How many alignments does Boyer-Moore try when matching the string 
    # GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu 
    # sequences) to the excerpt of human chromosome 1?  
    # (Don't consider reverse complements.)
    p_bm = BoyerMoore(PATTERN_1)
    occurrences_bm, num_aligments_bm, num_character_comparisons_bm =\
        hw.boyer_moore_with_counts(PATTERN_1, p_bm, TEXT)
    
    print(f"The naive exact matching algorithm tries {num_aligments_bm } alignments")
    
    # Q4
    # How many times does the string "GGCGCGGTGGCTCACGCCTGTAAT",  which is 
    # derived from a human Alu sequence, occur with up to 2 substitutions in the 
    # excerpt of human chromosome 1?  (Don't consider reverse complements here.)

    # Q5
    # How many total index hits are there when searching for occurences with up
    # to 2 substitutions in the excerpt of human chromosome 1?
    
    
    PATTERN_2 = "GGCGCGGTGGCTCACGCCTGTAAT"
    n = 2
    k = 8
    
    matches, hits = hw.index_approximate_matching(PATTERN_2, TEXT, 2, 8)
    print(f"The {len(matches)} matches found where {matches} "
          f"for a total of {hits} hits")
    
    # Q6
    # Write a function that, given a length-24 pattern P and given a SubseqIndex
    # object built with k = 8 and ival = 3, finds all approximate occurrences of 
    # P within T with up to 2 mismatches.
    # When using this function, how many total index hits are there when searching 
    # for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of 
    # human chromosome 1?  (Again, don't consider reverse complements.)
    
    n = 2
    k = 8
    interval = 3
    
    occurrences, n_hits = hw.subseq_approximate_matching(PATTERN_2, TEXT, n, k, interval)
    
    print(f"The number of total hits is {n_hits}")