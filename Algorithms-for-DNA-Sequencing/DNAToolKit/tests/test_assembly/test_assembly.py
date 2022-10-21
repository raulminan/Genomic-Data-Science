import os
import pytest

import time
from DNAToolKit.DNAToolKit import assembly, read_files

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
    
    def test_overlap_map_keys(self):
        reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
        map, _ = assembly.overlap_map(reads, 3)

        assert list(map.keys()) == [
            ('ABCDEFG', 'EFGHIJ'), 
            ('EFGHIJ', 'HIJABC'), 
            ('HIJABC', 'ABCDEFG')]

    def test_overlap(self):
        reads = ["ACGGATGATC", "GATCAAGT", "TTCACGGA"]
        map_, _ = assembly.overlap_map(reads, 3)

        assert map_ == {
            ("ACGGATGATC", "GATCAAGT"): 4,
            ("TTCACGGA", "ACGGATGATC"): 5,
        }
    
    def test_overlap_is_fast(self):
        reads, _ = read_files.read_fastq(
            os.path.join(
                os.path.dirname(__file__),
                "data",
                "ERR266411_1.for_asm.fastq"
            )
        )

        start = time.perf_counter()
        map, _ = assembly.overlap_map(reads, 30)
        end = time.perf_counter()
        time_ = end - start

        assert time_ < 20 # function shouldn't take "much more than 15s"

    def test_scs(self):
        set = ["ACGGTACGAGC", "GAGCTTCGGA", "GACACGG"]
        scs, _ = assembly.scs(set)
        
        assert scs == "GACACGGTACGAGCTTCGGA"
    
    def test_scs_list(self):
        set = ['ABC', 'BCA', 'CAB']
        _, scs_list = assembly.scs(set)

        assert scs_list == ['ABCAB', 'BCABC', 'CABCA']
    
    def test_greedy_scs(self):
        set1 = ["ABC", "BCA", "CAB"]
        set2 = ["ABCD", "CDBC", "BCDA"]

        assert assembly.greedy_scs(set1, 2) == "CABCA"
        assert assembly.greedy_scs(set2, 1) == "CDBCABCDA"
    
    def test_de_bruijn(self):
        read = "ACGCGTCG"
        k = 3

        nodes, edges = assembly.de_bruijn(read, k)

        assert nodes == {"GC", "AC", "GT", "CG", "TC"}
        assert edges == [
            ("AC", "CG"),
            ("CG", "GC"),
            ("GC", "CG"),
            ("CG", "GT"),
            ("GT", "TC"),
            ("TC", "CG"),
        ]
