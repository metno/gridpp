from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_barnes(self):
        structure = gridpp.BarnesStructure(1000)
        p1 = gridpp.Point(0, 0, 0, 0, gridpp.Cartesian)
        p2 = gridpp.Point(1000, 0, 0, 0, gridpp.Cartesian)
        self.assertEqual(1, structure.corr(p1, p1))
        self.assertEqual(0.6065306663513184, structure.corr(p1, p2))


if __name__ == '__main__':
    unittest.main()
