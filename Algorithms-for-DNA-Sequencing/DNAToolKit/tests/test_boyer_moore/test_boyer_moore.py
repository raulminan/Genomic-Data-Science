import pytest
import os

from DNAToolKit.DNAToolKit import boyer_moore, BoyerMoore

class Tests:
    """tests for boyer_moore module"""
    def test_bm_exact_matching(self):
        pattern = "word"
        text = "there would have been a time for such a word"
        lowercase_alphabet = "abcdefghijklmnopqrstuvwxyz "
        bm = BoyerMoore(pattern, lowercase_alphabet)
        occurences = boyer_moore.bm_exact_matching(pattern, text, bm)
        
        print(occurences)
        assert occurences == [40]

    def test_bm_exact_matching_with_counts(self):
        pattern = 'needle'
        text = 'needle need noodle needle'
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz ' 
        bm = BoyerMoore(pattern, lowercase_alphabet)
        
        occurrences, n_alignments, n_char_comparisons = boyer_moore.bm_exact_matching_with_counts(
            pattern, text, bm
        )

        assert occurrences == [0, 19]
        assert n_alignments == 5
        assert n_char_comparisons == 18

    def test_bm_approximate_matching(self):
        text = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
        pattern = 'to-morrow and to-morrow '
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz- ' 
        max_mismatches = 2

        occurrences = boyer_moore.bm_approximate_matching(
            pattern,
            text,
            max_mismatches,
            alphabet=lowercase_alphabet
        )
    
        assert occurrences == [0, 14]