from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_invalid_h(self):
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        for structure in structures:
            for h in [-1, np.nan]:
                with self.assertRaises(Exception) as e:
                    s = structure(h)
                with self.assertRaises(Exception) as e:
                    s = structure(h, 100)

    def test_barnes_h(self):
        hs = [-1, np.nan]
        for h in hs:
            with self.subTest(h=h):
                with self.assertRaises(Exception) as e:
                    structure = gridpp.BarnesStructure(h)

    def test_invalid_elevation(self):
        """Check that point elevations are ignored if one is missing"""
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        h = 2000
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 0, 0, gridpp.Cartesian)
        p3 = gridpp.Point(1000, 0, float('nan'), 0, gridpp.Cartesian)
        for structure in structures:
            with self.subTest(structure=structure):
                s1 = structure(h, 0)
                s2 = structure(h, 100)
                self.assertAlmostEqual(s1.corr(p1, p3), s1.corr(p1, p2))
                self.assertAlmostEqual(s2.corr(p1, p3), s2.corr(p1, p2))

    def test_invalid_v(self):
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        h = 2000
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 0, 0, gridpp.Cartesian)
        for structure in structures:
            for v in [-1, np.nan]:
                with self.assertRaises(Exception) as e:
                    s = structure(h, v)

    def test_invalid_w(self):
        structures = [gridpp.BarnesStructure, gridpp.CressmanStructure]
        h = 2000
        v = 100
        for structure in structures:
            for w in [-1, np.nan]:
                with self.assertRaises(Exception) as e:
                    s = structure(h, v, w)

    def test_invalid_cv(self):
        barnes = gridpp.BarnesStructure(2000)
        for dist in [-1, np.nan]:
            with self.assertRaises(Exception) as e:
                structure = gridpp.CrossValidation(barnes, dist)

    def test_multiple_structure(self):
        s1 = gridpp.CressmanStructure(2000, 2000, 2000)
        s2 = gridpp.CressmanStructure(200, 200, 200)
        s3 = gridpp.CressmanStructure(2, 2, 2)
        structure = gridpp.MultipleStructure(s1, s2, s3)

        expected = 0.6 # np.exp(-1)

        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 0, 0, gridpp.Cartesian)
        corr = structure.corr(p1, p2)
        self.assertAlmostEqual(corr, expected)

        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(0, 0, 100, 0, gridpp.Cartesian)
        corr = structure.corr(p1, p2)
        self.assertAlmostEqual(corr, expected)

        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(0, 0, 0, 1, gridpp.Cartesian)
        corr = structure.corr(p1, p2)
        self.assertAlmostEqual(corr, expected)

        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 100, 1, gridpp.Cartesian)
        corr = structure.corr(p1, p2)
        self.assertAlmostEqual(corr, expected**3)

    def test_multiple_structure_corr_vec(self):
        s1 = gridpp.CressmanStructure(5000, 11, 22)
        s2 = gridpp.CressmanStructure(33, 200, 44)
        s3 = gridpp.CressmanStructure(55, 66, 2)
        structure = gridpp.MultipleStructure(s1, s2, s3)

        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(0, 2500, 0, 0, gridpp.Cartesian)
        p3 = gridpp.Point(0, 2500, 100, 1, gridpp.Cartesian)
        corr = structure.corr(p1, p2)
        self.assertAlmostEqual(corr, 0.6)
        corr = structure.corr(p1, p3)
        self.assertAlmostEqual(corr, 0.6**3)

        # The only way to test the vector version of corr in python, is to go through OI. This is
        # because std::vector<Point> isn't exposed in python.
        N = 3
        y = [0, 0, 0]
        x = [0, 0, 0]
        z = [0, 0, 100]
        laf = [0, 0, 1]
        grid = gridpp.Points(y, x, z, laf, gridpp.Cartesian)
        points = gridpp.Points([0], [2500], [0], [0], gridpp.Cartesian)
        pratios = [1]
        pobs = [1]
        background = np.zeros([N])
        pbackground = [0]
        max_points = 10
        output = gridpp.optimal_interpolation(grid, background, points, pobs, pratios, pbackground, structure, max_points)
        np.testing.assert_array_almost_equal(output, [0.3, 0.3, 0.6**3/2])

    def test_clone(self):
        s1 = gridpp.CressmanStructure(5000, 11, 22)
        s2 = gridpp.CressmanStructure(33, 200, 44)
        s3 = gridpp.CressmanStructure(55, 66, 2)
        structure = gridpp.MultipleStructure(s1, s2, s3)
        structure_clone = structure.clone()
        del structure

        # Check that clone still works
        p1 = gridpp.Point(0, 0)
        p2 = gridpp.Point(0, 0)
        structure_clone.corr(p1, p2)

    def test_same_corr_after_clone(self):
        h = 850
        v = 92
        w = 0.44
        structures = [gridpp.BarnesStructure(h, v, w), gridpp.CressmanStructure(h, v, w)]
        structures += [gridpp.MultipleStructure(gridpp.BarnesStructure(1.3*h, v, w),
            gridpp.BarnesStructure(h, 1.3*v, w), gridpp.BarnesStructure(h, v, 1.3*w))]
        structures += [gridpp.CrossValidation(gridpp.BarnesStructure(h, v, w), 1000)]
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(500, 0, 50, 0.25, gridpp.Cartesian)
        for structure in structures:
            with self.subTest(structure=structure):
                structure_clone = structure.clone()
                ans = structure.corr(p1, p2)
                ans_clone = structure_clone.corr(p1, p2)
                self.assertEqual(ans, ans_clone)

                ans = structure.corr_background(p1, p2)
                ans_clone = structure_clone.corr_background(p1, p2)
                self.assertEqual(ans, ans_clone)


if __name__ == '__main__':
    unittest.main()
