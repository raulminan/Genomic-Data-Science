from Index import Index, SubseqIndex
from bm_preproc import BoyerMoore

class HomeWorkWeek2:
    """Contains functions to complete Week 2 homework"""
    def read_genome(filename):
        genome = ''
        with open(filename, 'r') as f:
            for line in f:
                # ignore header line with genome information
                if not line[0] == '>':
                    genome += line.rstrip()
        return genome

    def boyer_moore_with_counts(p, p_bm, t):
        """Do Boyer-Moore matching

        Args:
            p (str): pattern to match
            p_bm (bm_preproc.BoyerMoore): BoyerMoore object
            t (str): text to match

        Returns:
            list: list of occurences
            int: number of alingments
            int: number of character comparisons
        """
        i = 0
        occurrences = []
        num_alignments = 0
        num_character_comparisons = 0
        while i < len(t) - len(p) + 1:
            shift = 1
            mismatched = False
            num_alignments += 1
            for j in range(len(p)-1, -1, -1):
                num_character_comparisons += 1
                if p[j] != t[i+j]:
                    skip_bc = p_bm.bad_character_rule(j, t[i+j])
                    skip_gs = p_bm.good_suffix_rule(j)
                    shift = max(shift, skip_bc, skip_gs)
                    mismatched = True
                    break
            if not mismatched:
                occurrences.append(i)
                skip_gs = p_bm.match_skip()
                shift = max(shift, skip_gs)
            i += shift
        return occurrences, num_alignments, num_character_comparisons


    def naive_with_counts(p, t):
        """Do naive exact matching matching

        Args:
            p (str): pattern to match
            t (str): text to match

        Returns:
            list: list of occurences
            int: number of alingments
            int: number of character comparisons
        """
        occurrences = []
        num_alignments = 0
        num_character_comparisons = 0
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            num_alignments += 1
            for j in range(len(p)):  # loop over characters
                num_character_comparisons += 1
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
        return occurrences, num_alignments, num_character_comparisons

    def naive_2mm(p, t):
        """Do the naive exact match algorithms allowing up to 2 mismatches

        Args:
            p (str): pattern to match
            t (str): text to match the pattern to

        Returns:
            list: list with the indeces of occurences of the pattern p in text t
        """
        occurrences = set() 
        
        for i in range(len(t) - len(p) + 1):
            mismatches = 0 # loop over alignments
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                occurrences.add(i)  # all chars matched; record
        return list(occurrences)

    def boyer_approximate__matching(p, t, n):
        """Do approximate marching using the Boyer-Moore algorith and the pigeon
        principle

        Args:
            p (str): pattern to match
            t (str): text to match the pattern to
            n (int): maximum amount of mismatches allowed
        
        """
        segment_length = round(len(p) // (n+1))
        all_matches = set()
    
        for i in range(n+1):
            start = i*segment_length
            end = min((i+1)*segment_length, len(p))
            p_bm = BoyerMoore(p[start:end], alphabet="ACGT")
            matches, _, _ = HomeWorkWeek2.boyer_moore_with_counts(p[start:end], p_bm, t)
            
            for m in matches:
                if m - start < 0 or m - start + len(p) > len(t):
                    continue
                
                mismatches = 0
                # verify to the left
                for j in range(0, start):
                    if p[j] != t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break  
                        
                # verify to the right
                for j in range(end, len(p)):
                    if p[j] != t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break
            
                if mismatches <= n:
                    all_matches.add(m - start)
        
        return list(all_matches), len(matches)
        
    def index_approximate_matching(p, t, n, k):
        """Do approximate matching using the pigeonhole principle

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
        """Do approximate matching using subsequences

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
