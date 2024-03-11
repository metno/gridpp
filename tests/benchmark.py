import time
import numpy as np
import argparse
import collections
import gridpp

def main():
    parser = argparse.ArgumentParser(description='Runs gridpp benchmarks for processing performance')
    parser.add_argument('-j', type=int, help='Do a scaling test, by running on multiple cores (>= 2)', dest='num_cores')
    parser.add_argument('-s', type=float, default=1, help='Enlarge the inputs by this scaling factor to run a bigger test', dest='scaling')
    parser.add_argument('-n', type=int, default=1, help='Number of iterations to average over', dest='iterations')
    parser.add_argument('-t', help='Run only this function', dest="functions", nargs='*')

    args = parser.parse_args()

    if args.num_cores is not None and args.num_cores < 2:
        raise Exception("Error: Number of cores must be 2 or more")

    ss_time = time.time()

    np.random.seed(1000)

    version_major = int(gridpp.version().split('.')[0])
    version_minor = int(gridpp.version().split('.')[1])

    # Commonly used inputs
    structure = gridpp.BarnesStructure(10000)
    structure_grid = Structure(Grid(100, args.scaling).compute(), 10000)
    radius = 7
    quantile = 0.5
    thresholds = np.linspace(0, 1, 11)
    x = Input([2000, 2000, 10], args.scaling, sort_axis=2)
    y = Input([2000, 2000, 10], args.scaling, sort_axis=2)
    I1000 = Input([1000, 1000], args.scaling)
    Z1000 = Input([1000, 1000], args.scaling, random=False)
    I2000 = Input([2000, 2000], args.scaling)
    G1000 = Grid(1000, args.scaling)
    cache = Cache()

    """
    Define benchmarks, including their inputs and their expected runtime. For the inputs, use
    Computable objects to prevent expensive objects to be generated for all tests. These are then
    computed on the fly. If any other objects are used, these will be generated when the benchmark
    is defined. The goal here is to reduce the time it takes to initialize this program.

    The key is (function name to be tested, extra information for the test)
    The value is a dictionary with the following keys:
        expected (float): Expected runtime
        args (list): List of input arguments to the test
    """
    run = collections.OrderedDict()
    run[("Grid", "1000²")] = { "expected": 0.74, "args": (I1000, I1000)}
    run[("neighbourhood", "10000²")] = {"expected": 2.05, "args":(np.zeros([10000, 10000]), radius, gridpp.Mean)}
    run[("neighbourhood","2000² max")] = {"expected": 0.99, "args":(I2000, radius, gridpp.Max)}
    run[("neighbourhood_quantile_fast", "2000²")] = {"expected": 1.23, "args":(I2000, quantile, radius, thresholds)}
    run[("neighbourhood_quantile", "500²")] = {"expected": 1.70, "args":(Input([500, 500], args.scaling), quantile, radius)}
    run[("bilinear", "1000²")] = {"expected": 1.68, "args":(G1000, G1000, I1000)}
    run[("bilinear", "1000² x 50")] = {"expected": 4.42, "args":(G1000, G1000, Input([50, 1000, 1000], args.scaling))}
    run[("nearest", "1000²")] = {"expected": 1.52, "args":(G1000, G1000, I1000)}
    run[("nearest", "1000² x 50")] = {"expected": 1.93, "args":(G1000, G1000, Input([50, 1000, 1000], args.scaling))}
    run[("gridding", "200² 100000")] = {"expected": 0.61, "args":(Grid(200, args.scaling), Points(100000, args.scaling), np.zeros([100000]), 5000, 1, gridpp.Mean)}
    run[("gridding_nearest", "200² 100000")] = {"expected": 0.11, "args":(Grid(200, args.scaling), Points(100000, args.scaling), np.zeros([100000]), 1, gridpp.Mean)}
    run[("optimal_interpolation", "100² 1000")] = {"expected": 0.80, "args":(Grid(100, args.scaling), Input([100, 100], args.scaling), Points(1000, args.scaling), np.zeros(1000), np.ones(1000), np.ones(1000), structure, 20)}
    if version_minor >= 8:
        # This is slow in versions before 0.8.0
        run[("optimal_interpolation", "var len scale")] = {"expected": 0.91, "args":(Grid(100, args.scaling), Input([100, 100], args.scaling), Points(1000, args.scaling), np.zeros(1000), np.ones(1000), np.ones(1000), structure_grid, 20)}
    run[("dewpoint", "1e7")] = {"expected": 0.53, "args":(np.zeros(10000000) + 273.15, np.zeros(10000000))}
    run[("fill", "1e5")] = {"expected": 1.96, "args":(Grid(200, args.scaling), np.zeros([200, 200]), Points(100000, args.scaling), np.ones(100000) * 5000, 1, False)}
    run[("doping_square", "1e5")] = {"expected": 0.12, "args":(Grid(200, args.scaling), np.zeros([200, 200]), Points(100000, args.scaling), np.ones(100000) * 1, np.ones(100000, 'int') * 5, False)}
    run[("doping_circle", "1e5")] = {"expected": 2.00, "args":(Grid(200, args.scaling), np.zeros([200, 200]), Points(100000, args.scaling), np.ones(100000) * 1, np.ones(100000) * 5000, False)}
    run[("local_distribution_correction", "")] = {"expected": 1.31, "args":(Grid(200, args.scaling), np.zeros([200, 200]), Points(1000, args.scaling), np.ones(1000) * 1, np.ones(1000) * 1, structure, 0.1, 0.9, 5)}
    run[("full_gradient", "1000²")] = {"expected": 1.59, "args": (G1000, G1000, I1000, I1000, I1000)}
    if version_major > 0 or version_minor > 6:
        run[("calc_gradient", "2000²")] = {"expected": 0.45, "args": (np.random.rand(2000, 2000) * 100, np.zeros([2000,2000]), gridpp.LinearRegression, 10, 0, 100, 0)}
        run[("mask_threshold_downscale_consensus", "")] = {"expected": 0.91, "args":(Grid(100, args.scaling), G1000, np.random.rand(100, 100, 10), np.random.rand(100, 100, 10), np.random.rand(1000, 1000, 10), np.random.rand(1000, 1000), gridpp.Lt, gridpp.Mean)}
    run[("neighbourhood_search", "2000² 7x7")] = {"expected": 1.11, "args": (np.random.rand(2000, 2000), np.random.rand(2000, 2000), 3, 0.7, 1, 0.1, np.random.rand(2000, 2000) < 0.5)}
    run[("window", "1000²")] = {"expected": 1.67, "args": (np.random.rand(100000, 1000), 101, gridpp.Mean, False, False)}
    run[("gamma_inv", "5*201*476")] = {"expected": 1.168, "args": (np.random.rand(5*201*476)*0.9 + 0.05, np.random.rand(5*201*476), np.random.rand(5*201*476))}
    run[("apply_curve", "")] = {"expected": 0.06, "args": (I2000, np.random.rand(2000), np.random.rand(2000), gridpp.OneToOne, gridpp.OneToOne)}
    run[("apply_curve", "gridded")] = {"expected": 0.87, "args": (I2000, x, y, gridpp.OneToOne, gridpp.OneToOne)}
    run[("test_vec3_input")] = {"expected": 0.35, "args": (np.zeros([2000, 2000, 10], np.float32),)}
    run[("init_vec3")] = {"expected": 0.8, "args": (1000, 1000, 200)}

    if args.num_cores is not None:
        print("Gridpp parallelization test (gridpp version %s)" % gridpp.version())
    else:
        print("Gridpp benchmark (gridpp version %s)" % gridpp.version())
        print("Expected results from Intel i7 3.40 Ghz")
    print("-----------------------------------------------------------------")
    if args.num_cores is not None:
        print("Function                               1 core %2d cores  Scaling" % args.num_cores)
    else:
        print("Function                             Expected     Time     Diff")
    num_cores = [1]
    if args.num_cores is not None:
        num_cores += [args.num_cores]

    for key in run.keys()       :
        try:
            timings = dict()
            for num_core in num_cores:
                timings[num_core] = 0

            if isinstance(key, tuple):
                name = key[0] + " " + str(key[1])
                func = eval("gridpp." + key[0])
            else:
                name = key
                func = eval("gridpp." + key)
            if args.functions is not None:
                if func.__name__ not in args.functions:
                    continue

            # Allow functions to fail (useful when benchmarking older versions of gridpp
            # where functions may not be defined).
            for num_core in num_cores:
                gridpp.set_omp_threads(num_core)
                for it in range(args.iterations):
                    # Evaluate Computable inputs
                    input_args = list()
                    for i, k in enumerate(run[key]["args"]):
                        if issubclass(type(k), Computable):
                            k = cache.add(k)
                        elif isinstance(k, list):
                            inner = list()
                            for j in k:
                                if issubclass(type(j), Computable):
                                    j = cache.add(j)
                                inner += k
                            k = inner
                        input_args += [k]

                    # Run the actual test
                    s_time = time.time()
                    func(*input_args)
                    e_time = time.time()

                    curr_time = e_time - s_time
                    timings[num_core] += curr_time
        except Exception as e:
            raise e
            print("Could not run", key, e)
            continue

        for num_core in num_cores:
            timings[num_core] /= args.iterations

        if args.num_cores is None:
            diff = (timings[1] - run[key]["expected"] * args.scaling) / (run[key]["expected"]  * args.scaling) * 100
            string = "%-36s %8.2f %8.2f %8.2f %%" % (name, run[key]["expected"] * args.scaling, timings[1], diff)
        else:
            maximum_speedup = args.num_cores
            minimum_speedup = 1
            actual_speedup = timings[1] / timings[args.num_cores]

            # Compute the fraction between minimum and maximum speedup
            scaling = ((actual_speedup - minimum_speedup) / (maximum_speedup - minimum_speedup))

            # A better way to compute scaling? Computed in time space, instead of in factor space
            # Probably not a fair way, since scaling will stay fixed for higher number of cores
            # expected_time = timings[1] / args.num_cores
            # scaling = 1 - (timings[args.num_cores] - expected_time) / (timings[1] - expected_time)

            string = "%-36s %8.2f %8.2f %8.2f %%" % (name, timings[1], timings[args.num_cores], scaling * 100)

            # factor = timings[1] / timings[args.num_cores]
            # string = "%-36s %8.2f %8.2f %8.2f" % (name, timings[1], timings[args.num_cores], factor)
        print(string)


def get_shape(obj):
    shape = None
    try:
        shape = obj.shape
    except Exception as e:
        pass
    try:
        shape = obj.size()
    except Exception as e:
        pass
    return shape


class Cache:
    """Stores data in a cache:
    cache[object type][hashable key]
    """
    def __init__(self):
        self.cache = dict()

    def has(self, obj):
        if type(obj) not in self.cache:
            return False
        s = obj in self.cache[type(obj)]
        return s

    def add(self, obj):
        """Adds object to cache and return the computed value of the object. If it already exists,
        the computed object is returned"""
        if not self.has(obj):
            if not issubclass(type(obj), Computable):
                raise ValueError("Only Computable objects can be added to cache")

            if type(obj) not in self.cache:
                self.cache[type(obj)] = dict()

            ret = obj.compute()
            self.cache[type(obj)][obj] = ret
        else:
            ret = self.cache[type(obj)][obj]

        return ret


class Computable:
    def compute(self):
        raise NotImplementedError()

    def __neq__(self, other):
        return self != other

class Input(Computable):
    def __init__(self, shape, scaling, sort_axis=None, random=False):
        self.shape = shape
        self.scaling = scaling
        self.sort_axis = sort_axis
        self.random = random

    def compute(self):
        shape = [int(self.shape[0] * self.scaling)] + self.shape[1:]
        if self.random:
            output = np.random.rand(*shape)*10
        else:
            output = np.zeros(shape)

        if self.sort_axis is not None:
            output = np.sort(output, axis=self.sort_axis)

        return output

    def __eq__(self, other):
        return self.shape == other.shape and self.scaling == other.scaling

    def __str__(self):
        return "Input " + ','.join([str(s) for s in self.shape])

    def __hash__(self):
        return hash(str(self.shape))

class Grid(Computable):
    def __init__(self, shape, scaling):
        self.shape = shape
        self.scaling = scaling

    def compute(self):
        return gridpp.Grid(*np.meshgrid(np.linspace(0, 1, self.shape), np.linspace(0, 1, int(self.shape * self.scaling))))

    def __eq__(self, other):
        return self.shape == other.shape and self.scaling == other.scaling

    def __str__(self):
        return "Grid " + f"{self.shape}"

    def __hash__(self):
        return self.shape

class Points(Computable):
    def __init__(self, shape, scaling):
        self.shape = shape
        self.scaling = scaling

    def compute(self):
        return gridpp.Points(np.random.rand(self.shape), np.random.rand(self.shape))

    def __eq__(self, other):
        return self.shape == other.shape and self.scaling == other.scaling

    def __str__(self):
        return f"Points {self.shape}"

    def __hash__(self):
        return self.shape

class Structure(Computable):
    def __init__(self, grid, h):
        self.grid = grid
        self.h = h

    def compute(self):
        q = np.ones(self.grid.size())
        return gridpp.BarnesStructure(self.grid, self.h*q, 0*q, 0*q)

    def __eq__(self, other):
        return self.h == other.h and self.grid.size() == other.grid.size()

    def __str__(self):
        return f"Grid {self.grid.size()}"

    def __hash__(self):
        return hash(str(self.grid.size()))

if __name__ == "__main__":
    main()
