from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_rotated(self):
        A = gridpp.Point(-8.59, -8.89)
        B = gridpp.Point(-3.41, -11.89)
        C = gridpp.Point(2.60, -1.5)
        D = gridpp.Point(-2.60, 1.5)
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(-3, -5.20)))

    def test_not_in_order(self):
        """Check case where 4 points are not in order"""
        A = gridpp.Point(0, 0)
        B = gridpp.Point(1, 1)
        C = gridpp.Point(0, 1)
        D = gridpp.Point(1, 0)
        self.assertFalse(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0.6, 0.6)))

    def test_skew(self):
        A = gridpp.Point(0, 0)
        B = gridpp.Point(1, 0.25)
        C = gridpp.Point(1, 1.25)
        D = gridpp.Point(0, 1)
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0.5, 0.5)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0, 0)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(1, 1.25)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0.5, 0.25)))
        self.assertFalse(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0.5, 0)))

    def test_cw(self):
        A = gridpp.Point(0, 0)
        B = gridpp.Point(1, 0)
        C = gridpp.Point(1, 1)
        D = gridpp.Point(0, 1)
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0.5, 0.5)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0, 0)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(1, 1)))
        self.assertFalse(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(-0.1, -0.1)))

    def test_ccw(self):
        A = gridpp.Point(0, 0)
        B = gridpp.Point(0, 1)
        C = gridpp.Point(1, 1)
        D = gridpp.Point(1, 0)
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0.5, 0.5)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(0, 0)))
        self.assertTrue(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(1, 1)))
        self.assertFalse(gridpp.point_in_rectangle(A, B, C, D, gridpp.Point(-0.1, -0.1)))


if __name__ == '__main__':
    unittest.main()
