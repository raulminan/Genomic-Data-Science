import pytest
from DNAToolKit import alignment

class Tests:
    def test_global_alignment_transition(self):
        x = "AGTCGATCGAGCTAGCGCA"
        y = "GGTCGATCGAGCTAGCGCA"
        
        assert alignment.global_alignment(x, y) == 2

    def test_global_alignment_tranversion(self):
        x = "AGTCGATCGAGCTAGCGCA"
        y = "TGTCGATCGAGCTAGCGCA"

        assert alignment.global_alignment(x, y) == 4

    def test_global_alignment_skip(self):
        x = "AGTCGATCGAGCTAGCGCA"
        y = "AGTCGATCGGCTAGCGCA"

        assert alignment.global_alignment(x, y) == 8

    def test_global_alignment_all(self):
        x = "AGTCGATCGAGCTAGCGCA"
        y = "TGTCGATCAGCTAGCGTA"

        assert alignment.global_alignment(x, y) == 14