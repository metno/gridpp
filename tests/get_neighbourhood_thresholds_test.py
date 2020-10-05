from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]

"""Simple check
20 21 22 23 24
15 16 17 18 19
10 11 12 13 nan
5  6  7  nan  9
0  1  2  3  4
"""
values = np.reshape(range(25), [5, 5]).astype(float)
values[1, 3] = np.nan
values[2, 4] = np.nan
values = np.array(values)

class Test(unittest.TestCase):
    def test_invalid_arguments(self):
        """Check that exception is thrown for invalid arguments"""
        field = np.ones([5, 5])
        nums = [-1, 0]
        for num in nums:
            with self.assertRaises(ValueError) as e:
                gridpp.get_neighbourhood_thresholds(field, num)

    def test_high_num(self):
        """Check proper results when num > size of the input"""
        field = np.reshape(np.arange(4), [2, 2])
        nums = [4, 5, 6]
        for num in nums:
            q = gridpp.get_neighbourhood_thresholds(field, num)
            np.testing.assert_array_equal(q, [0, 1, 2, 3])

    def test_duplicates(self):
        """Check a real example with duplicates"""
        field = np.reshape([0, 0, 2, 3, 4, 5, 6, 11, 8, 9, 10, 11], [3, 4])
        nums = range(1, 5)
        for num in nums:
            q = gridpp.get_neighbourhood_thresholds(field, num)
            self.assertTrue((q >= 0).all())
            self.assertTrue((q <= 11).all())


if __name__ == '__main__':
    unittest.main()
