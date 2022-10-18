from typing import Tuple
from .BoyerMoore import BoyerMoore

def bm_exact_matching(pattern: str, text: str, boyer_moore: BoyerMoore) -> list:
    """Does Boyer-Moore matching

    Parameters
    ----------
    pattern : str
        pattern to match against text
    text : str
        text to match pattern to 
    boyer_moore : BoyerMoore
        BoyerMoore object

    Returns
    -------
    list
        list with all the offsets of text where pattern occurs
    """
    i = 0
    occurrences = []
    while i < len(text) - len(pattern) + 1:
        shift = 1
        mismatched = False
        for j in range(len(pattern)-1, -1, -1):
            if pattern[j] != text[i+j]:
                skip_bc = boyer_moore.bad_character_rule(j, text[i+j])
                skip_gs = boyer_moore.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = boyer_moore.match_skip()
            shift = max(shift, skip_gs)
        i += shift

    return occurrences

def bm_exact_matching_with_counts(pattern: str,
                                  text: str,
                                  boyer_moore: BoyerMoore) -> Tuple[list, int, int]:
    """Does Boyer-Moore matching and counts the number of alignments and 
    characters compared

    Parameters
    ----------
    pattern : str
        pattern to match against text
    text : str
        text to match pattern to 
    boyer_moore : BoyerMoore
        BoyerMoore object

    Returns
    -------
    Tuple[list, int, int]
        list: list of offsets where pattern occurs in text
        int: number of alignments
        int: number of character comparisons
    """
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(text) - len(pattern) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(pattern)-1, -1, -1):
            num_character_comparisons += 1
            if pattern[j] != text[i+j]:
                # skip_bc = BoyerMoore.bad_character_rule(j, text[i+j])
                # skip_gs = BoyerMoore.good_suffix_rule(j)
                skip_bc = boyer_moore.bad_character_rule(j, text[i+j])
                skip_gs = boyer_moore.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            # skip_gs = BoyerMoore.match_skip()
            skip_gs = boyer_moore.match_skip()
            shift = max(shift, skip_gs)
        i += shift

    return occurrences, num_alignments, num_character_comparisons

def bm_approximate_matching(pattern: str, 
                            text: str, 
                            max_mismatches: int,
                            alphabet: str = "ATCG") -> list:
    """Does approximate Boyer-Moore matching using the pigeon hole principle

    Parameters
    ----------
    pattern : str
        pattern to match
    text : str
        text to match the pattern to
    max_mismatches : int
        maximum amount of mismatches allowed
    alphabet: str
        alphabet to build BoyerMoore objects, by default "ATCG"

    Returns
    -------
    list
        list with all offsets where pattern occurs in text
    """
    segment_length = round(len(pattern) // (max_mismatches+1))
    all_matches = set()

    for i in range(max_mismatches+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(pattern))
        bm = BoyerMoore(pattern[start:end], alphabet=alphabet)
        matches = bm_exact_matching(pattern[start:end], text, bm)
        
        for m in matches:
            if m - start < 0 or m - start + len(pattern) > len(text): 
                # check that it doesn't run off the beggining or past the end
                continue
            
            mismatches = 0
            # verify to the left
            for j in range(0, start):
                if pattern[j] != text[m-start+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break  
                    
            # verify to the right
            for j in range(end, len(pattern)):
                if pattern[j] != text[m-start+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
        
            if mismatches <= max_mismatches:
                all_matches.add(m - start)

    return list(all_matches)