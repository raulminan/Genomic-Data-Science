import os
import pytest

from DNAToolKit.DNAToolKit import naive

class Tests:
    "tests for the naive module"
    def test_naive_matching(self):
        pattern = "CTGT"
        ten_as = "AAAAAAAAAA"
        text = ten_as + "CTGT" + ten_as + "CTTT" + ten_as + "CGGG" + ten_as
        mismatches = 2

        occurrences = naive.naive_matching(pattern, text, mismatches)

        assert occurrences == [10, 24, 38]
