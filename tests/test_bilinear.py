from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]


class BilinearTest(unittest.TestCase):
    def test_neighbourhood(self):
        """Simple check
        50  6  7  8
        40  3  4  5
        30  0  1  2
            0 10 20
        """
        lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)
        values[:] = np.reshape(range(9), lons.shape)
        points = gridpp.Points([25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55], [-1, -1, -1, 0, 0, 0, 10, 10, 10, 20, 20, 20, 21, 21, 21])
        output = gridpp.nearest(grid, points, values)
        np.testing.assert_array_equal(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))
        output = gridpp.bilinear(grid, points, values)
        np.testing.assert_array_equal(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))

    def test_incompatible_sizes(self):
        """ Check for exception when input is different size than grid"""
        lats = [[0, 0, 0, 0], [1, 1, 1, 1]]
        lons = [[0, 1, 2, 3], [0, 1, 2, 3]]
        grid = gridpp.Grid(lats, lons)  # 2 x 4 grid
        points = gridpp.Points([0, 1], [0, 1])
        for N in [1, 3, 5]:
            values = np.zeros([N, N])
            with self.assertRaises(Exception) as e:
                gridpp.bilinear(grid, grid, values)
            with self.assertRaises(Exception) as e:
                gridpp.bilinear(grid, points, values)

    def test_rotation(self):
        """Check that the interpolation gives the same result when grid is rotated"""
        values = np.array([[0, 0, 0, 0], [0, 0, 1, 0], [0, 2, 3, 0], [0, 0, 0, 0]])
        rotations = range(0, 360, 10)
        lons0 = np.array([[-24, -12, 0, 12], [-24, -12, 0, 12], [-24, -12, 0, 12], [-24, -12, 0, 12]])
        lats0 = np.array([[-6, -6, -6, -6], [-3, -3, -3, -3], [3, 3, 3, 3], [6, 6, 6, 6]])
        lons_p0 = -6
        lats_p0 = 0
        lons = np.zeros([4, 4])
        lats = np.zeros([4, 4])
        for rotation in rotations:
            angle = rotation * 2 * 3.14159265 / 360
            lons_p = lons_p0 * np.cos(angle) - lats_p0 * np.sin(angle)
            lats_p = lons_p0 * np.sin(angle) + lats_p0 * np.cos(angle)
            points = gridpp.Points([lats_p], [lons_p])
            # print("Rotation: %d" % rotation)

            for i in range(lons0.shape[0]):
                for j in range(lons0.shape[1]):
                    lons[i, j] = lons0[i, j] * np.cos(angle) - lats0[i, j] * np.sin(angle)
                    lats[i, j] = lons0[i, j] * np.sin(angle) + lats0[i, j] * np.cos(angle)
            if 0:
                import matplotlib.pylab as mpl
                mpl.plot(lons, lats, '-ko')
                mpl.plot(lons_p, lats_p, '-ro')
                mpl.show()
            # mpl.plot(lons, lats)
            # mpl.show()
            # sys.exit()
            grid = gridpp.Grid(lats, lons)
            output = gridpp.bilinear(grid, points, values)[0]
            self.assertAlmostEqual(output, 1.5, 5)

    def test_parallelogram(self):
        """Skew a square"""
        values = np.array([[0, 1], [2, 3]])
        for skew in [-1,-0.1, 0, 0.1,1]:
            """Move the bottom edge from left to right"""
            lats = [[-1, -1], [1, 1]]
            lons = [[-1+skew, 1+skew], [-1-skew, 1-skew]]
            grid = gridpp.Grid(lats, lons)
            output = gridpp.bilinear(grid, gridpp.Points([0], [0]), values)[0]
            self.assertEqual(output, 1.5)

            """Move the left edge from bottom to top"""
            lats = [[-1+skew, -1-skew], [1+skew, 1-skew]]
            lons = [[-1, 1], [-1, 1]]
            grid = gridpp.Grid(lats, lons)
            output = gridpp.bilinear(grid, gridpp.Points([0], [0]), values)[0]
            self.assertEqual(output, 1.5)

    def test_non_parallelogram(self):
        """Non-parallelogram quadrilateral"""
        values = np.array([[0, 1], [2, 3]])
        lats = [[0, 0], [1, 2]]
        lons = [[0, 1], [0.1, 1]]
        grid = gridpp.Grid(lats, lons)
        output = gridpp.bilinear(grid, gridpp.Points([0.5], [0.5]), values)[0]
        # TODO: This value has not been confirmed by an independent calculation
        self.assertAlmostEqual(output, 1.1570623, 5)

    # Check that this doesn't fail
    def test_vertical_parallel(self):
        values = np.array([[0, 1], [2, 3]])
        lons = np.transpose([[-10, -9], [-11, -10]])
        lats = np.transpose([[54, 54], [55, 55.1]])
        lons = [[-10, -11], [-9, -10]]
        lats = [[54, 55], [54, 55.1]]
        grid = gridpp.Grid(lats, lons)
        output = gridpp.bilinear(grid, gridpp.Points([55], [-10]), values)[0]

        # Trygveasp's Ireland example
        values = np.array([[0, 1], [2, 3]])
        lons = np.transpose([[-10.2027884, -9.97075176], [-10.2363253, -10.0035467]])
        lats = np.transpose([[54.8651619, 54.884182], [54.9986267, 55.0177078]])
        grid = gridpp.Grid(lats, lons)
        output = gridpp.bilinear(grid, gridpp.Points([55], [-10]), values)[0]

    def test_with_missing_values(self):
        """Check that nearest neighbour is returned when a missing value is in the box"""
        values = np.array([[0, np.nan, 2, 3], [0, 1, 2, 3]])
        lats = [[0, 0, 0, 0], [1, 1, 1, 1]]
        lons = [[0, 1, 2, 3], [0, 1, 2, 3]]
        grid = gridpp.Grid(lats, lons)
        # The first three points have a missing value in its box, whereas the 4th does not
        output = gridpp.bilinear(grid, gridpp.Points([0, 0.3, 0.9, 0.5], [0, 0.3, 0.9, 2.5]), values)
        np.testing.assert_array_equal([0, 0, 1, 2.5], output)

    def test_grid_to_grid(self):
        """Check grid to grid itnerpolation"""
        values = np.reshape(np.arange(4), [2, 2])
        lons1, lats1 = np.meshgrid([0, 1], [0, 1])
        lons2, lats2 = np.meshgrid([0, 0.5, 1], [0, 0.5, 1])
        grid1 = gridpp.Grid(lats1, lons1)
        grid2 = gridpp.Grid(lats2, lons2)
        output = gridpp.bilinear(grid1, grid2, values)
        np.testing.assert_array_equal([[0, 0.5, 1], [1, 1.5, 2], [2, 2.5, 3]], output)

    def test_grid_to_grid_3d(self):
        """Check grid to grid itnerpolation"""
        T = 3
        values = np.reshape(np.arange(4), [2, 2])
        values = np.expand_dims(values, 0)
        values = np.repeat(values, T, axis=0)
        lons1, lats1 = np.meshgrid([0, 1], [0, 1])
        lons2, lats2 = np.meshgrid([0, 0.5, 1], [0, 0.5, 1])
        grid1 = gridpp.Grid(lats1, lons1)
        grid2 = gridpp.Grid(lats2, lons2)
        output = gridpp.bilinear(grid1, grid2, values)
        for t in range(T):
            np.testing.assert_array_equal([[0, 0.5, 1], [1, 1.5, 2], [2, 2.5, 3]], output[t, ...])

    def test_dimensions_mismatch(self):
        lons, lats = np.meshgrid([0, 1], [0, 1])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros([3, 1, 1])
        with self.assertRaises(Exception) as e:
            gridpp.bilinear(grid, grid, values)

    def test_empty_output_points(self):
        lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)

        points = gridpp.Points([], [])
        output = gridpp.bilinear(grid, points, values)
        np.testing.assert_array_equal(output, [])

    def test_empty_input_grid(self):
        igrid = gridpp.Grid([[]], [[]])
        points = gridpp.Points([0, 5, 10], [0, 5, 10])
        output = gridpp.bilinear(igrid, points, np.zeros([0, 0]))
        np.testing.assert_array_equal(output, [np.nan, np.nan, np.nan])

    def test_empty_output_grid(self):
        lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
        igrid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)

        grid = gridpp.Grid([[]], [[]])
        output = gridpp.bilinear(igrid, grid, values)
        np.testing.assert_array_equal(output, [[]])

    def test_weird(self):
        x = np.reshape([-117.3010559, -116.8340607, -117.2735291, -116.803299], [2, 2]).transpose()
        y = np.reshape([57.71018982, 57.69472504, 57.95975876, 57.94418716], [2, 2]).transpose()
        x0 = -117.299057 - x[0][0]
        y0 = 57.71475601 - y[0][0]
        x = x - x[0][0]
        y = y - y[0][0]
        values = np.reshape(np.arange(4), [2, 2]).transpose()

        for do_transpose in [False, True]:
            with self.subTest(structure=do_transpose):
                x = x.transpose()
                y = y.transpose()
                values = values.transpose()
                grid = gridpp.Grid(y, x)
                points = gridpp.Points([y0], [x0])
                q = gridpp.bilinear(grid, points, values)
                self.assertAlmostEqual(q[0], 0.04017778, 6)

    def check(self):
        N = 4
        lons, lats = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))
        lons1, lats1 = np.meshgrid(np.linspace(0, 1, 101), np.linspace(0, 1, 101))
        grid = gridpp.Grid(lats, lons)
        grid1 = gridpp.Grid(lats1, lons1)
        values = np.random.rand(lons.shape[0], lons.shape[1])
        print(lons.shape, values.shape)
        values1 = gridpp.bilinear(grid, grid1, values)
        import matplotlib.pylab as mpl
        mpl.pcolormesh(lons1, lats1, values1)
        mpl.show()



if __name__ == '__main__':
    unittest.main()
