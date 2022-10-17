"""functions for matching using index data structures"""
from .Index import Index
from .SubseqIndex import SubseqIndex

def substring_approximate_matching(p, t, n, k):
    """Do approximate matching using substrings and the pigeonhole principle

    Args:
        p (str): pattern to match
        t (str): text to match the pattern to
        n (int): maximum amount of mismatches allowed
        k (int): k-mer to use for indexing
    """
    segment_length = round(len(p) // (n+1))
    all_matches = set()
    idx = Index(t, k)
    n_hits = 0

    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        hits = idx.query(p[start:end])
        n_hits += len(hits)
        
        for hit in hits:
            if hit - start < 0 or hit - start + len(p) > len(t):
                continue
            
            mismatches = 0
            # verify to the left
            for j in range(0, start):
                if not p[j] == t[hit-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break  
                    
            # verify to the right
            for j in range(end, len(p)):
                if not p[j] == t[hit-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
        
            if mismatches <= n:
                all_matches.add(hit - start)
    
    return list(all_matches), n_hits

def subseq_approximate_matching(p, t, n, k, ival):
    """Do approximate matching using subsequences and the pigeonhole principle

    Args:
        p (str): pattern to match
        t (str): text to match the pattern to
        n (int): maximum amount of mismatches allowed
        k (int): k-mer to use for indexing
        ival (int): interval to use for create sub sequences
    """
    all_matches = set()
    idx = SubseqIndex(t, k, ival)
    n_hits = 0
    
    for i in range(n+1):
        start = i 
        hits = idx.query(p[start:])
        n_hits += len(hits)
        for hit in hits:
            if hit - start < 0 or hit - start + len(p) > len(t):
                continue
            
            mismatches = 0
            t
            for j in range(0, len(p)):
                if not p[j] == t[hit-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break  
        
            if mismatches <= n:
                all_matches.add(hit - start)
    
    return list(all_matches), n_hits
    