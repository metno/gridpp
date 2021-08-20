from __future__ import print_function
import unittest
import gridpp
import numpy as np
import os


class Test(unittest.TestCase):
    def test_missing_values(self):
        self.assertTrue(np.isnan(gridpp.sea_level_pressure(np.nan, 20, 290)))
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, 101325, np.nan, 290)
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, 101325, 20, np.nan)

    def test_unphysical_values(self):
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, -1, 20, 290)
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, 101325, 20, -1)
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, 101325, 20, 290, -1)
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, 101325, 20, 290, 2)
        self.assertRaises(RuntimeError, gridpp.sea_level_pressure, 101325, 20, 290, 0.7, -1)

    def test_vector(self):
        ps = [101315., 101000., 102300., 99513.]
        alt = [38, 34, 51, 69]
        temperature = [290, 273, 293, 295]
        rh = [0.1, 0.5, 0.8, 0.9]
        dewpoint = [np.nan, np.nan, np.nan, np.nan]
        truth = [gridpp.sea_level_pressure(ps[i], alt[i], temperature[i], rh[i], dewpoint[i]) for i in range(len(ps))]
        np.testing.assert_array_almost_equal(truth, gridpp.sea_level_pressure(ps, alt, temperature, rh, dewpoint))

    def test_dimension_mismatch(self):
        """Check for exception when vector arguments have different sizes"""
        ps = [101315., 101000., 102300., 99513.]
        alt = [38, 34, 51, 69]
        temperature = [290, 273, 293, 295]
        rh = [0.1, 0.5, 0.8, 0.9]
        dewpoint = [np.nan, np.nan, np.nan, np.nan]
        with self.assertRaises(Exception) as e:
            gridpp.sea_level_pressure(ps + [100], alt, temperature, rh, dewpoint)
        with self.assertRaises(Exception) as e:
            gridpp.sea_level_pressure(ps, alt + [100], temperature, rh, dewpoint)
        with self.assertRaises(Exception) as e:
            gridpp.sea_level_pressure(ps, alt, temperature + [100], rh, dewpoint)
        with self.assertRaises(Exception) as e:
            gridpp.sea_level_pressure(ps, alt, temperature, rh + [100], dewpoint)
        with self.assertRaises(Exception) as e:
            gridpp.sea_level_pressure(ps, alt, temperature, rh, dewpoint +[100] )

    def test_pressure_results(self):
            self.assertAlmostEqual(gridpp.sea_level_pressure(101325.0, 20, 273.15), 101578.0)
            self.assertAlmostEqual(gridpp.sea_level_pressure(101325.0, 50, 273.15), 101960.25)
  

if __name__ == '__main__':
    unittest.main()
