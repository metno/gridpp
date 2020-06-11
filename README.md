# Gridded post-processor

![alt text](https://travis-ci.org/metno/gridpp.svg?branch=master "https://travis-ci.org/metno/gridpp") ![alt text](https://coveralls.io/repos/metno/gridpp/badge.svg?branch=master&service=github "https://coveralls.io/github/metno/gridpp?branch=master")

Gridpp a is post-processing tool for gridded weather forecasts. It consists of a **library** of commonly-used methods and a **command-line tool** that applies these methods to forecast fields in NetCDF files.

Gridpp is written in C++ but offers python bindings to the functions in the library. The tool is used at MET Norway to produce operational weather forecasts for Yr (https://www.yr.no).

Gridpp is currently under active development and the current version is a prototype for testing. Feedback
is welcome, either by using the issue tracker in Github, or by contacting Thomas Nipen (thomasn@met.no).

## Documentation
For information on how to use gridpp, check out the wiki at https://github.com/metno/gridpp/wiki.

## Features
- Methods for **downscaling** a forecast from a coarse grid to a fine grid
- Methods for **calibrating** a downscaled grid, such as quantile mapping
- Computationally efficient **neighbourhood** methods to compute neighbourhood min, mean, max, and any quantile.
- Data assimilation using **optimal interpolation** (OI) to merge observations and gridded forecasts (deterministic or ensemble)
- Efficient data structures for nearest location lookup in a vector or grid of locations
- Command-line client with support for Netcdf files with flexibility in how variables and dimensions are configured

## Example

The following computes a moving neighbourhood mean with a half-width of 7 gridpoints (i.e. 15x15 neighbourhood)

```python
import gridpp
import numpy as np

field = np.random.rand(300, 200)
halfwidth = 7
gridpp.neighbourhood(field, halfwidth, gridpp.Mean)
```

## Required dependencies
- [Boost](https://www.boost.org/) >= 1.59
- [Armadillo](http://arma.sourceforge.net/) >= 6.6
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [Netcdf](https://www.unidata.ucar.edu/software/netcdf/)

On Ubuntu Bionic, these can be installed like this:
```bash
sudo apt-get update
sudo apt-get install libboost-all-dev
sudo apt-get install libgsl0-dev libblas-dev
sudo apt-get install netcdf-bin libnetcdf-dev
sudo apt-get install cmake
sudo apt-get install libarmadillo6 libarmadillo-dev
```

Note that Ubuntu Xenial only has Armadillo 6.5 in its apt repository. In that case you need to install  [Armadillo 6.6](http://arma.sourceforge.net/) or later manually.

## Installing the python interface to gridpp only
The python package is available on [pypi.org](https://pypi.org/project/gridpp/). Provided you have installed the dependencies listed above, you can install the most recent release of the python package as follows:
```bash
pip3 install gridpp --user
```

To check that the installation worked, run the following in python3:
```python
import gridpp
print(gridpp.version())
```

## Full gridpp installation from source

1. Either download the source code from the [latest release](https://github.com/metno/gridpp/releases), unzip
   the file and navigate into the extracted folder; or clone the repo from github.

2. Create a build directory

```bash
mkdir build
```

3. Run cmake to set up installation

```bash
cd build
cmake ..
```

This will set up the gridpp command-line tool to be installed in /usr/local/bin on Ubuntu. To specify a custom installation path, use:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=<custom path>
```

4. Run cmake to install

```bash
cmake --build .
cmake --build . --target install
```

## Copyright and license
Copyright © 2014-2020 Norwegian Meteorological Institute. Gridpp is licensed under the GNU LEsser General
Public License (LGPL). See LICENSE file.

## Contact
E-mail: Thomas Nipen (thomasn@met.no)
