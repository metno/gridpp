#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.command.build import build as build_orig
from distutils.core import setup, Extension
from setuptools import setup, Extension
import glob
import itertools

__version__ = '${PROJECT_VERSION}'

# We need to swap the order of build_py and build_ext
# https://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module
from distutils.command.build import build
from setuptools.command.install import install
import distutils.command.install as orig

class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        # self.do_egg_install()
        orig.install.run(self)


module = Extension('_gridpp',
        sources=glob.glob('src/api/*.cpp') + glob.glob('src/api/*.c') + ['src/gridppPYTHON_wrap.cxx'],
        #language="c++",
        #swig_opts=['-I./include', '-c++', '-I/usr/include/python3.6m'],
        libraries=["gsl", "gslcblas", "proj", "armadillo"],
        library_dirs=["/usr/lib/x86_64-linux-gnu/"],
        include_dirs=['./include']
)

setup (
    name='gridpp',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='Post-processing software for gridded weather forecasts',
    long_description="Gridded post-processor\n======================\nGridpp is a command-line tool that post-processes weather forecasts in NetCDF format.  The program performs two types of post-processing: Downscaling and calibration. Gridpp downscales forecast from a coars grid to a finer grid using a variety of interpolation methods. Gridpp then calibrates the forecasts by applying corrections to each gridpoint. Gridpp is modular, so any combination of downscaling and calibration can be selected.\n\nFor information on how to use the software, check out the wiki page: https://github.com/metno/gridpp/wiki",

    # The project's main homepage.
    url='https://github.com/metno/gridpp',

    # Author details
    author='Thomas Nipen',
    author_email='thomasn@met.no',

    # Choose your license
    license='LGPL-3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Information Analysis',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    # What does your project relate to?
    keywords='meteorology quality control observation weather',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', '*tests*']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy>=1.7', 'scipy', 'six', 'future'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
    #    'dev': ['check-manifest'],
        'test': ['coverage', 'pep8'],
    #    'test': ['pytest'],
    },

    test_suite="gridpp.tests",
    ext_modules = [module],
    py_modules = ["gridpp"],
    # cmdclass={'build': build},
    cmdclass={'build': CustomBuild, 'install': CustomInstall},

    include_package_data=True,

)
