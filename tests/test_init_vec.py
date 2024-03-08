from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_vec2(self):
        X = -1
        for func in [gridpp.init_vec2, gridpp.init_ivec2]:
            np.testing.assert_array_equal(func(2, 2, X), [[X, X], [X, X]])
            np.testing.assert_array_equal(func(1, 2, X), [[X, X]])
            np.testing.assert_array_equal(func(2, 1, X), [[X], [X]])

    def test_vec3(self):
        X = -1
        for func in [gridpp.init_vec3, gridpp.init_ivec3]:
            np.testing.assert_array_equal(func(2, 2, 3, X), [[[X,X,X], [X,X,X]], [[X, X, X], [X,X,X]]])
            np.testing.assert_array_equal(func(1, 2, 3, X), [[[X,X,X], [X,X,X]]])
            np.testing.assert_array_equal(func(2, 1, 3, X), [[[X,X,X]], [[X, X, X]]])
            np.testing.assert_array_equal(func(2, 2, 1, X), [[[X], [X]], [[X], [X]]])

if __name__ == '__main__':
    unittest.main()
