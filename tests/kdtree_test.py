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
        tree = gridpp.KDTree([0, 1000, 2000], [0, 1000, 2000], gridpp.Cartesian)
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

    def test_rad2deg(self):
        self.assertAlmostEqual(gridpp.KDTree_rad2deg(1), 180 / 3.14159265, 5)
        self.assertAlmostEqual(gridpp.KDTree_rad2deg(-1), -180 / 3.14159265, 5)
        self.assertAlmostEqual(gridpp.KDTree_rad2deg(0), 0, 5)

    def test_calc_distance(self):
        config = list()
        # lat0, lon0, lat1, lon1, delta, expected
        config += [[60, 10, 61, 11, 0.1, 124080.79]]
        config += [[60, 10, 60, 10, 0, 0]]
        config += [[90, 10, -90, 10, 0, 20037508]]
        config += [[0, 0, 0, 180, 0, 20037508]]
        config += [[60.5, 5.25, -84.75, -101.75, 0, 16879114]]

        for c in config:
            self.assertEqual(len(c), 6)
            p0 = gridpp.Point(c[0], c[1])
            p1 = gridpp.Point(c[2], c[3])
            delta = c[4]
            expected = c[5]

            # Point version
            with self.subTest(config=c, type="Point"):
                self.assertAlmostEqual(expected, gridpp.KDTree.calc_distance(p0, p1), delta=delta)

            # Scalar version
            with self.subTest(config=c, type="Scalar"):
                self.assertAlmostEqual(expected, gridpp.KDTree.calc_distance(c[0], c[1], c[2], c[3]), delta=delta)

    def test_calc_straight_distance(self):
        config = list()
        # lat0, lon0, lat1, lon1, delta, expected
        config += [[60, 10, 61, 11, 2, 124080.79]]
        config += [[60, 10, 60, 10, 0, 0]]
        config += [[90, 10, -90, 10, 0, 6.378137e6*2]]
        config += [[0, 0, 0, 180, 10, 6.378137e6*2]]
        config += [[60.5, 5.25, -84.75, -101.75, 0, 12367265.0]]

        for c in config:
            self.assertEqual(len(c), 6)
            p0 = gridpp.Point(c[0], c[1])
            p1 = gridpp.Point(c[2], c[3])
            delta = c[4]
            expected = c[5]

            # Point version
            with self.subTest(config=c, type="Point"):
                self.assertAlmostEqual(expected, gridpp.KDTree.calc_straight_distance(p0, p1), delta=delta)

            # Scalar version
            with self.subTest(config=c, type="Scalar"):
                self.assertAlmostEqual(expected, gridpp.KDTree.calc_straight_distance(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z), delta=delta)

    def test_calc_distance_limit(self):
        p0 = gridpp.Point(0, 0)
        p1 = gridpp.Point(0.001, 0.001)
        self.assertAlmostEqual(157.42953491210938, gridpp.KDTree_calc_distance(0,0,0.001,0.001));
        self.assertAlmostEqual(157.42953491210938, gridpp.KDTree_calc_straight_distance(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z));

    def test_calc_distance_fast(self):
        config = list()
        #          lat0,lon0,lat1,lon1, delta, dist
        config += [[60,10, 60,10,  10, 0]]
        config += [[90,10, -90,10, 10, 20037508]]
        config += [[0,0,   0,180,  10, 20037508]]
        config += [[60,10, 61,11,  400, 124080.79]]
        config += [[89,0, 90,0,  10, 111319.62]]
        config += [[89,0, 90,180,  10, 111319.62]]
        config += [[89,0, 89.9,180,  6000, 111319.62]]

        for c in config:
            self.assertEqual(len(c), 6)
            delta = c[4]
            expected = c[5]

            self.assertAlmostEqual(expected, gridpp.KDTree.calc_distance_fast(c[0], c[1], c[2], c[3]), delta=delta)

    def test_calc_distance_fast_across_date_line(self):
        config = list()
        for lat in [-90, -89, 0, 89, 90]:
            config += [[lat, 180, lat, -180, 10, 0]]
        config += [[0, 179, 0, -179, 100, 222639.64]]
        for lat in [-90, -89, 0, 89]:
            config += [[lat, 180, lat+1, -180, 10, 111319.4921875]]

        for c in config:
            self.assertEqual(len(c), 6)
            p0 = gridpp.Point(c[0], c[1])
            p1 = gridpp.Point(c[2], c[3])
            delta = c[4]
            expected = c[5]

            self.assertAlmostEqual(expected, gridpp.KDTree.calc_distance_fast(c[0], c[1], c[2], c[3]), delta=delta)

    def test_radius_match(self):
        """Check that points right on the radius edge count as a match"""
        points = gridpp.Points([0, 1000, 2000], [0, 0, 0], [0, 0, 0], [0, 0, 0], gridpp.Cartesian)
        I = points.get_neighbours(900, 0, 501)
        np.testing.assert_array_equal(I, [1])
        I = points.get_neighbours(900, 0, 99.99)
        np.testing.assert_array_equal(I, [])
        I = points.get_neighbours(0, 0, 1000)
        np.testing.assert_array_equal(I, [0])
        I = points.get_neighbours(0, 0, 1001)
        np.testing.assert_array_equal(I, [0, 1])
        I = points.get_neighbours(0, 0, 1001, False)
        np.testing.assert_array_equal(I, [1])

    def test_empty_constructor(self):
        tree = gridpp.KDTree()
        self.assertEqual(tree.get_coordinate_type(), gridpp.Geodetic)

    def test_invalid_coords(self):
        lats = [91, -91, np.nan, 0]
        lons = [0, 0, 0, np.nan]
        for i in range(len(lats)):
            curr_lats = [lats[i]]
            curr_lons = [lons[i]]
            with self.subTest(lat=lats[i], lon=lons[i]):
                with self.assertRaises(ValueError) as e:
                    tree = gridpp.KDTree(curr_lats, curr_lons, gridpp.Geodetic)

    def test_valid_coords(self):
        lats = [90.000001, -90.0000001]
        lons = [0, 0]
        for i in range(len(lats)):
            curr_lats = [0, lats[i]]
            curr_lons = [0, lons[i]]
            with self.subTest(lat=lats[i], lon=lons[i]):
                tree = gridpp.KDTree(curr_lats, curr_lons, gridpp.Geodetic)
                I = tree.get_nearest_neighbour(0, 0)
                self.assertEqual(I, 0)

    def test_wrap_lon(self):
        """Check that longitudes outside [-180, 180] are correctly handled"""
        lons = [-360, 0, 360]
        for i in range(len(lons)):
            curr_lats = [0]
            curr_lons = [lons[i]]
            with self.subTest(lon=lons[i]):
                tree = gridpp.KDTree(curr_lats, curr_lons, gridpp.Geodetic)
                I, dist = tree.get_neighbours_with_distance(0, 0, 1e9)
                self.assertEqual(I[0], 0)
                self.assertAlmostEqual(dist[0], 0)

                I, dist = tree.get_neighbours_with_distance(0, 180, 1e9)
                self.assertEqual(I[0], 0)
                diameter_of_earth = 12756274.0
                self.assertAlmostEqual(dist[0], diameter_of_earth)


if __name__ == '__main__':
    unittest.main()
