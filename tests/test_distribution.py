from __future__ import print_function
import unittest
import gridpp
import numpy as np
import time
import collections


class Test(unittest.TestCase):
    def test_1(self):
        levels = [0.5, 0.5, 0.5]
        shape = [1, 2, 7.5]
        scale = [2, 2, 1]
        output = gridpp.gamma_inv(levels, shape, scale)
        expected = [1.386,3.357, 7.169]
        np.testing.assert_array_almost_equal(output, expected, 3)

    def test_invalid_levels(self):
        ok_args = collections.OrderedDict({
                'levels': [0.1],
                'shape': [1],
                'scale': [1],
        })
        invalid_args = {
                'levels': [[-0.1], [1.1], [np.nan]],
                'shape': [[-1], [np.nan]],
                'scale': [[-1], [np.nan]],
                }

        for key in invalid_args.keys():
            for arg in invalid_args[key]:
                args0 = ok_args.copy()
                args0[key] = arg
                q = [args0[f] for f in args0]
                with self.subTest(key=key, arg=arg):
                    with self.assertRaises(ValueError) as e:
                        output = gridpp.gamma_inv(*q)

        invalid = [-0.1, 1.1, np.nan]
        for i in range(len(invalid)):
            levels = [0.5, invalid[i], 0.5]
            shape = [1,1,1]
            scale = [1,1,1]
            with self.assertRaises(ValueError) as e:
                output = gridpp.gamma_inv(levels, shape, scale)



if __name__ == '__main__':
    unittest.main()
