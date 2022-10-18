"""functions for naive matching"""
from typing import Tuple

def naive_matching(pattern: str, text: str, max_mismatches: int = 0) -> list:

    """Performs naive matching allowing a given amount of mismatches

    Parameters
    ----------
    pattern : str
        pattern to match against text
    text : str
        text to match the pattern to
    max_mismatches : int, optional
        number of allowed mismatches, by default 0

    Returns
    -------
    list
        offsets where pattern matches text
    """
    occurrences = set()
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            occurrences.add(i)
        
    return sorted(list(occurrences))

def naive_matching_with_counts(
    pattern: str,
    text: str,
    max_mismatches: int = 0) -> Tuple[list, int, int]:
    """Performs naive matching allowing a given amount of mismatches and also
    returns the number of aligments and the number of character comparisons done

    Parameters
    ----------
    pattern : str
        pattern to match against text
    text : str
        text to match the pattern to
    max_mismatches : int, optional
        number of allowed mismatches, by default 0

    Returns
    -------
    Tuple[list, int, int]
        list: all offsets where pattern matches against text
        int: number of alignments
        int: number of char comparisons
    """
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(text) - len(pattern) + 1):  
        mismatches = 0
        num_alignments += 1
        for j in range(len(pattern)):  
            num_character_comparisons += 1
            if text[i+j] != pattern[j]:  
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            occurrences.append(i)  

    return occurrences, num_alignments, num_character_comparisons

# TODO implement naive_with_rc