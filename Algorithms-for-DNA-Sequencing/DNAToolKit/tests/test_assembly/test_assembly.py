import os
import pytest
from DNAToolKit import assembly

class Tests:
    """tests for assembly module"""
    def test_overlap_without(self):
        a = "TTACGT"
        b = "GTACCGT"

        assert assembly.overlap(a, b) == 0

    def test_overlap_with(self):
        a = "TTACGT"
        b = "CGTACCGT"

        assert assembly.overlap(a, b) == 3

    def test_naive_overlap_map(self):
        reads = ["ACGGATGATC", "GATCAAGT", "TTCACGGA"]
        map_ = assembly.naive_overlap_map(reads, 3)

        assert map_ == {
            ("ACGGATGATC", "GATCAAGT"): 4,
            ("TTCACGGA", "ACGGATGATC"): 5,
        }