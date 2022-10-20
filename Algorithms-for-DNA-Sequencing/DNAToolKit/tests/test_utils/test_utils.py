import pytest

from DNAToolKit.DNAToolKit import utils

class Tests:
    def test_edit_distance(self):
        x = "shake spea"
        y = "Shakespear"
        distance = utils.edit_distance(x, y)
        assert distance == 3
    
    def test_edit_distance(self):
        pattern = "GCGTATGC"
        text = "TATTGGCTATACGGTT"
        distance = utils.edit_distance_fewest_edits(pattern, text)

        assert distance == 2