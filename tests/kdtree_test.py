from __future__ import print_function
import unittest
import gridpp
import numpy as np


class KDTreeTest(unittest.TestCase):
    def test(self):
        lats = [60, 61, 62]
        lons = [10, 10, 12]
        tree = gridpp.KDTree(lats, lons)
        np.testing.assert_array_equal(tree.get_neighbours(60, 10, 1), [0])
        np.testing.assert_array_equal(tree.get_neighbours(60, 10, 112000), [0, 1])

if __name__ == '__main__':
    unittest.main()
