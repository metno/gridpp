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

    def test_flat(self):
        tree = gridpp.KDTree([0, 1000, 2000], [0, 1000, 2000], True)
        q, dist = tree.get_neighbours_with_distance(100, 100, 1000)
        self.assertEqual(dist, 100 * np.sqrt(2))

    """ Check that duplicate points are found """
    def test_duplicate_points(self):
        tree = gridpp.KDTree([50, 50, 51], [0, 0, 10])
        I = tree.get_neighbours(50, 0.001, 1000)
        self.assertEqual(len(I), 2)
        self.assertTrue(0 in I)
        self.assertTrue(1 in I)

    """ Check that duplicate points that are equal to lookup point are found """
    def test_duplicate_points_identical(self):
        tree = gridpp.KDTree([50, 50, 51], [0, 0, 10])
        I = tree.get_neighbours(50, 0, 1000)
        self.assertEqual(len(I), 2)
        self.assertTrue(0 in I)
        self.assertTrue(1 in I)
    def test_pole(self):
        tree = gridpp.KDTree([89, 89, 90], [0, 180, 0])
        I, dist = tree.get_neighbours_with_distance(90, 0, 1000)
        self.assertEqual(len(I), 1)
        self.assertAlmostEqual(dist[0], 0)
        I, dist = tree.get_neighbours_with_distance(90, 10, 1000)
        self.assertEqual(len(I), 1)
        self.assertAlmostEqual(dist[0], 0)
        I, dist = tree.get_neighbours_with_distance(89.999, 0, 1000)
        self.assertEqual(len(I), 1)

    def test_pole_with_duplicate_points(self):
        tree = gridpp.KDTree([89, 89, 90, 90], [0, 180, 0, 10])
        I, dist = tree.get_neighbours_with_distance(90, 0, 1000)
        self.assertEqual(len(I), 2)
        self.assertTrue(2 in I)
        self.assertTrue(3 in I)
        self.assertAlmostEqual(dist[0], 0)
        self.assertAlmostEqual(dist[1], 0)

if __name__ == '__main__':
    unittest.main()
