from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_basic(self):
        x = [0, 1000, 2000, 3000, np.nan]
        structure_corr = dict()
        barnes = gridpp.BarnesStructure(2000)
        structure_corr[barnes] = [1, 0.8824968934059143, 0.6065306663513184, 0.32465246319770813, 0]
        structure_corr[gridpp.CressmanStructure(2000)] = [1, 0.6, 0, 0, 0]
        structure_corr[gridpp.CrossValidation(barnes, 1000)] = [0, 0, 0.6065306663513184, 0.32465246319770813, 0]
        N = len(x)
        for structure, corr in structure_corr.items():
            is_cv = isinstance(structure, gridpp.CrossValidation)
            for i in range(N):
                with self.subTest(structure=type(structure), i=i):
                    p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
                    p2 = gridpp.Point(x[i], 0, 0, 0, gridpp.Cartesian)
                    funcs = [structure.corr, structure.corr_background]
                    if is_cv:
                        funcs = [structure.corr_background]
                    for func in funcs:
                        self.assertAlmostEqual(corr[i], func(p1, p2))

                        # Check that changing the order does not change the results
                        self.assertAlmostEqual(corr[i], func(p2, p1))

                        # Check identical points
                        if not is_cv and not np.isnan(x[i]):
                            self.assertAlmostEqual(1, func(p2, p2))

    def test_invalid_h(self):
        hs = [-1, np.nan]
        for h in hs:
            with self.subTest(h=h):
                with self.assertRaises(Exception) as e:
                    structure = gridpp.BarnesStructure(h)

    def test_invalid_hmax(self):
        with self.assertRaises(Exception) as e:
            structure = gridpp.BarnesStructure(2000, 100, 0, -1)

    def test_spatial(self):
        y = [[0, 0]]
        x = [[0, 2500]]
        grid = gridpp.Grid(y, x, y, y, gridpp.Cartesian)
        h = [[2500, 1]]
        v = [[0, 0]]
        l = [[0, 0]]
        min_rho = 0.1
        structure = gridpp.BarnesStructure(grid, h, v, l, min_rho)
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(0, 2500, 0, 0, gridpp.Cartesian)

        # Length scale at p1 is 2500
        expected = np.sqrt(-2*np.log(min_rho)) * 2500
        self.assertAlmostEqual(expected, structure.localization_distance(p1), 4)
        self.assertAlmostEqual(0.6, structure.corr(p1, p2), 1)

        # Length scale at p2 is 1
        expected = np.sqrt(-2*np.log(min_rho)) * 1
        self.assertAlmostEqual(expected, structure.localization_distance(p2), 4)
        self.assertAlmostEqual(0, structure.corr(p2, p1), 1)

    def test_spatial_invalid_arguments(self):
        """Check dimension mismatch"""
        y, x = np.meshgrid(np.linspace(0, 1, 2), np.linspace(0,1,3))
        grid = gridpp.Grid(y, x, y, y, gridpp.Cartesian)

        valid = np.ones([3, 2])
        invalid = [np.ones([3,4]), np.ones([2, 2]), np.ones([2, 4])]
        structure = gridpp.BarnesStructure(grid, valid, valid, valid)

        for inval in invalid:
            with self.assertRaises(ValueError) as e:
                structure = gridpp.BarnesStructure(grid, inval, valid, valid)
            with self.assertRaises(ValueError) as e:
                structure = gridpp.BarnesStructure(grid, valid, inval, valid)
            with self.assertRaises(ValueError) as e:
                structure = gridpp.BarnesStructure(grid, valid, valid, inval)

    def test_hmax(self):
        hmaxs = [0, 1000, 2000, 10000]
        p0 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        dist_ans = {0:1, 1000:0.8824968934059143, 2000:0.6065306663513184, 3000:0.32465246319770813}
        for hmax in hmaxs:
            for dist, ans in dist_ans.items():
                with self.subTest(hmax=hmax, dist=dist):
                    structure = gridpp.BarnesStructure(2000, 0, 0, hmax)
                    corr = structure.corr(p0, gridpp.Point(dist, 0, 0, 0, gridpp.Cartesian))
                    if dist > hmax:
                        self.assertEqual(0, corr)
                    else:
                        self.assertEqual(ans, corr)


if __name__ == '__main__':
    unittest.main()
