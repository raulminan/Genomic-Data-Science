import pytest

from DNAToolKit import utils

class Tests:
    def test_edit_distance(self):
        x = "shake spea"
        y = "Shakespear"
        distance = utils.edit_distance(x, y)
        assert distance == 3