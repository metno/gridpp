from __future__ import print_function
import unittest
import gridpp
import numpy as np
import psutil
import os


def memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


class SwigTest(unittest.TestCase):
    def test_int_output(self):
        """Checks if ints are passed back
        """
        output = gridpp.test_ivec_output()
        self.assertEqual(type(output[0]), np.int32)
        output = gridpp.test_ivec2_output()
        self.assertEqual(type(output[0][0]), np.int32)

    def test_float_output(self):
        output = gridpp.test_vec_output()
        self.assertEqual(type(output[0]), np.float32)
        output = gridpp.test_vec2_output()
        self.assertEqual(type(output[0][0]), np.float32)
        output = gridpp.test_vec3_output()
        self.assertEqual(type(output[0][0][0]), np.float32)


if __name__ == '__main__':
    unittest.main()
