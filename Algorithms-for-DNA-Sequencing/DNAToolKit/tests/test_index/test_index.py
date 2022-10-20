import os
import pytest

from DNAToolKit.DNAToolKit import index, read_files, naive

class Tests:
    """Tests for index module"""

    @pytest.fixture
    def setup(self):
        self.file = os.path.join(os.path.dir(__file__), "data", "chr1.GRCh38.excerpt.fasta")

    def test_substring_approximate_matching_no_mismatches(self):
        text = "GCTACGATCTAGAATCTA"
        pattern = "TCTA"
        max_mismatches = 0
        kmer_length = 4

        occurences, hits = index.substring_approximate_matching(
            pattern, 
            text, 
            max_mismatches, 
            kmer_length
        )

        occurences_naive = naive.naive_matching(pattern, text, 0)
        assert occurences == occurences_naive
    
    def test_substring_approximate_matching_2_mismatches(self):
        text = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
        pattern = 'to-morrow and to-morrow '
        max_mismatches = 2
        kmer_length = 8

        occurrences, _ = index.substring_approximate_matching(
            pattern, text, max_mismatches, kmer_length
        )
        
        occurrences_naive = naive.naive_matching(
            pattern, text, max_mismatches
        )

        assert occurrences == occurrences_naive


    def test_subsequence_approximate_matching(self):
        text = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
        pattern = 'to-morrow and to-morrow '
        max_mismatches = 2
        kmer_length = 8
        interval = 3

        occurences, n_hits = index.subseq_approximate_matching(
            pattern,
            text,
            max_mismatches,
            kmer_length,
            interval
        )

        assert occurences == [0, 14]
        assert n_hits == 6
 