"""functions for matching using index data structures"""
from .Index import Index
from .SubseqIndex import SubseqIndex



def substring_approximate_matching(pattern: str,
                                   text: str,
                                   max_mismatches: int,
                                   kmer_length: int) -> list[int]:
    """Does approximate matching building an Index data structure allowing a 
    given number of mismatches

    Parameters
    ----------
    pattern : str
        pattern to match
    text : str
        text to match pattern to
    max_mismatches : int
        max number of mismatches allowed
    kmer_length : int
        kmer_length to use to build the Index data scructure

    Returns
    -------
    list[int]
        list with the offsets in text where pattern occurs
    """
    segment_length = len(pattern) // (max_mismatches+1)
    all_matches = set()
    idx = Index(text, kmer_length)
    n_hits = 0
    
    for i in range(max_mismatches + 1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(pattern))
        hits = idx.query(pattern[start:end])
        n_hits += len(hits)
        
        for hit in hits:
            if (hit - start < 0) or (hit - start + len(pattern) > len(text)):
                continue
            
            mismatches = 0
            # verify to the left
            for j in range(0, start):
                if pattern[j] != text[hit-start+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break  
                    
            # verify to the right
            for j in range(end, len(pattern)):
                if pattern[j] != text[hit-start+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
        
            if mismatches <= max_mismatches:
                all_matches.add(hit - start)
    
    return sorted(list(all_matches)), n_hits

# TODO fix this function
def subseq_approximate_matching(pattern, text, max_mismatches, kmer_length, interval):
    """Do approximate matching using subsequences
    Args:
        p (str): pattern to match
        t (str): text to match the pattern to
        n (int): maximum amount of mismatches allowed
        k (int): k-mer to use for indexing
        ival (int): interval to use for create sub sequences
    """
    segment_length = len(pattern) // (max_mismatches+1)
    all_matches = set()
    idx = SubseqIndex(text, kmer_length, interval)
    print(idx.index)
    n_hits = 0
    
    for i in range(max_mismatches + 1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(pattern))
        hits = idx.query(pattern[start:end])
        n_hits += len(hits)
        
        for hit in hits:
            if (hit - start < 0) or (hit - start + len(pattern) > len(text)):
                continue
            
            mismatches = 0
            # verify to the left
            for j in range(0, start):
                if pattern[j] != text[hit-start+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break  
                    
            # verify to the right
            for j in range(end, len(pattern)):
                if pattern[j] != text[hit-start+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
        
            if mismatches <= max_mismatches:
                all_matches.add(hit - start)
    
    return sorted(list(all_matches)), n_hits
