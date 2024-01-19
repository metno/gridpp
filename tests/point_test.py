from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_initialize(self):
        lat = 60
        lon = 0
        elev = 10
        laf = 0.3
        p = gridpp.Point(lat, lon, elev, laf, gridpp.Geodetic)
        self.assertAlmostEqual(p.lat, lat)
        self.assertAlmostEqual(p.lon, lon)
        self.assertAlmostEqual(p.elev, elev)
        self.assertAlmostEqual(p.laf, laf, delta=0.00001)

        s, x, y, z = gridpp.convert_coordinates(lat, lon, gridpp.Geodetic)
        self.assertAlmostEqual(p.x, x)
        self.assertAlmostEqual(p.y, y)
        self.assertAlmostEqual(p.z, z)

        p = gridpp.Point(lat, lon, elev, laf, gridpp.Geodetic, x, y, z)
        self.assertAlmostEqual(p.lat, lat)
        self.assertAlmostEqual(p.lon, lon)
        self.assertAlmostEqual(p.elev, elev)
        self.assertAlmostEqual(p.laf, laf, delta=0.00001)
        self.assertAlmostEqual(p.x, x)
        self.assertAlmostEqual(p.y, y)
        self.assertAlmostEqual(p.z, z)


if __name__ == '__main__':
    unittest.main()
