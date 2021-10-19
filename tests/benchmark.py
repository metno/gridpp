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
    parser.add_argument('-t', help='Run only this function', dest="function")

    args = parser.parse_args()

    if args.num_cores is not None and args.num_cores < 2:
        raise Exception("Error: Number of cores must be 2 or more")

    input = dict()
    grids = dict()
    points = dict()
    np.random.seed(1000)

    for i in [10, 50, 100, 200, 500, 1000, 2000, 10000]:
        input[i] = np.random.rand(int(i * args.scaling), i)*10
    for i in [10, 50, 100, 200, 500, 1000]:
        grids[i] = gridpp.Grid(*np.meshgrid(np.linspace(0, 1, i), np.linspace(0, 1, int(i * args.scaling))))
    for i in [1000, 100000]:
        # points[i] = gridpp.Points(np.linspace(0, 1, i), np.zeros(i))
        points[i] = gridpp.Points(np.random.rand(i) * 10, np.random.rand(i) * 10)
    structure = gridpp.BarnesStructure(10000)
    radius = 7
    quantile = 0.5
    thresholds = np.linspace(0, 1, 11)
    run = collections.OrderedDict()
    run[("Grid", "1000²")] = {"expected": 0.74, "args":np.meshgrid(np.linspace(0, 1, int(1000 * args.scaling)), np.linspace(0, 1, 1000))}
    run[("neighbourhood", "10000²")] = {"expected": 2.05, "args":(np.zeros([10000, 10000]), radius, gridpp.Mean)}
    run[("neighbourhood","2000² max")] = {"expected": 0.99, "args":(input[2000], radius, gridpp.Max)}
    run[("neighbourhood_quantile_fast", "2000²")] = {"expected": 1.23, "args":(input[2000], quantile, radius, thresholds)}
    run[("neighbourhood_quantile", "500²")] = {"expected": 1.70, "args":(input[500], quantile, radius)}
    run[("bilinear", "1000²")] = {"expected": 1.68, "args":(grids[1000], grids[1000], input[1000])}
    run[("bilinear", "1000² x 50")] = {"expected": 4.42, "args":(grids[1000], grids[1000], np.repeat(np.expand_dims(input[1000], 0), 50, axis=0))}
    run[("nearest", "1000²")] = {"expected": 1.52, "args":(grids[1000], grids[1000], input[1000])}
    run[("nearest", "1000² x 50")] = {"expected": 2.30, "args":(grids[1000], grids[1000], np.repeat(np.expand_dims(input[1000], 0), 50, axis=0))}
    run[("gridding", "1000² 100000")] = {"expected": 0.53, "args":(grids[1000], points[100000], np.zeros([100000]), 5000, 1, gridpp.Mean)}
    run[("gridding_nearest", "1000² 100000")] = {"expected": 0.13, "args":(grids[1000], points[100000], np.zeros([100000]), 1, gridpp.Mean)}
    run[("optimal_interpolation", "1000² 1000")] = {"expected": 1.57, "args":(grids[1000],
        input[1000], points[1000], np.zeros(1000), np.ones(1000), np.ones(1000), structure, 20)}
    run[("dewpoint", "1e7")] = {"expected": 0.53, "args":(np.zeros(10000000) + 273.15, np.zeros(10000000))}
    run[("fill", "1e5")] = {"expected": 0.52, "args":(grids[1000], np.zeros([1000, 1000]),
        points[100000], np.ones(100000) * 5000, 1, False)}
    run[("doping_square", "1e5")] = {"expected": 0.16, "args":(grids[1000], np.zeros([1000, 1000]),
        points[100000], np.ones(100000) * 1, np.ones(100000, 'int') * 5, False)}
    run[("doping_circle", "1e5")] = {"expected": 0.52, "args":(grids[1000], np.zeros([1000, 1000]),
        points[100000], np.ones(100000) * 1, np.ones(100000) * 5000, False)}
    run[("local_distribution_correction", "")] = {"expected": 0.52, "args":(grids[1000], np.zeros([1000, 1000]),
        points[1000], np.ones(1000) * 1, np.ones(1000) * 1, structure, 0.1, 0.9, 5)}
    run[("full_gradient", "1000²")] = {"expected": 1.59, "args": (grids[1000], grids[1000], np.zeros([1000,1000]), np.zeros([1000,1000]), np.zeros([1000,1000]))}
    run[("calc_gradient", "1000²")] = {"expected": 0.18, "args": (np.random.rand(1000, 1000) *
        100, np.zeros([1000,1000]), gridpp.LinearRegression, 10, 0, 0, 0)}
    run[("window", "1000²")] = {"expected": 0.78, "args": (np.random.rand(1000, 1000), 100, gridpp.Mean, False, False)}

    if args.num_cores is not None:
        print("Gridd parallelization test (gridpp version %s)" % gridpp.version())
    else:
        print("Gridd benchmark (gridpp version %s)" % gridpp.version())
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
            if args.function is not None:
                if func.__name__ != args.function:
                    continue

            # Allow functions to fail (useful when benchmarking older versions of gridpp
            # where functions may not be defined).
            for num_core in num_cores:
                gridpp.set_omp_threads(num_core)
                for it in range(args.iterations):
                    s_time = time.time()
                    func(*run[key]["args"])
                    e_time = time.time()
                    curr_time = e_time - s_time
                    timings[num_core] += curr_time
        except Exception as e:
            print("Could not run", key, e)
            continue

        for num_core in num_cores:
            timings[num_core] /= args.iterations

        if args.num_cores is None:
            diff = (timings[1] - run[key]["expected"] * args.scaling) / (run[key]["expected"]  * args.scaling) * 100
            string = "%-36s %8.2f %8.2f %8.2f %%" % (name, run[key]["expected"] * args.scaling, timings[1], diff)
        else:
            scaling = timings[1] / timings[args.num_cores] / args.num_cores
            expected = timings[1] / args.num_cores
            scaling = 1 - (timings[args.num_cores] - expected) / (timings[1] - expected)
            # scaling = (1 - timings[args.num_cores] / timings[1]) * (args.num_cores + 1)

            string = "%-36s %8.2f %8.2f %8.2f %%" % (name, timings[1], timings[args.num_cores], scaling * 100)
        print(string)


if __name__ == "__main__":
    main()
