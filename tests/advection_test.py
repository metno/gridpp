from __future__ import print_function
import unittest
import gridpp
import numpy as np


class AdvectionTest(unittest.TestCase):
    def test_1(self):
        Y = 3
        X = 3
        y_dist = np.ones([Y, X])
        for x in range(3):
            y_dist[:, x] = x
        x_dist = np.zeros([Y, X])
        y, x = gridpp.advection_implicit(y_dist, x_dist, 1)


if __name__ == '__main__':
    unittest.main()
