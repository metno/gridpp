from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_get_box(self):
        # 2x3 grid
        lats = [[0, 0, 0], [1, 1, 1]]
        lons = [[0, 1, 2], [0.25, 1.25, 2.25]]
        grid = gridpp.Grid(lats, lons)
        # Check that the point lat=0.4,lon=1.25 is in between lat[0:1] and lon[1:2]
        np.testing.assert_array_equal(grid.get_box(0.4, 1.25), [True, 0, 1, 1, 2])
        np.testing.assert_array_equal(grid.get_box(0.4, 1), [True, 0, 0, 1, 1])


if __name__ == '__main__':
    unittest.main()
