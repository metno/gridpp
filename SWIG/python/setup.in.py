#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import setuptools.command.install
import shutil
from distutils.sysconfig import get_python_lib
from distutils.command.build import build as build_orig
from distutils.core import Extension
import itertools

__version__ = '${PROJECT_VERSION}'

class CompiledLibInstall(setuptools.command.install.install):
    """
    Specialized install to install to python libs
    """

    def run(self):
        """
        Run method called by setup
        :return:
        """
        # Get filenames from CMake variable
        filenames = '${PYTHON_INSTALL_FILES}'.split(';')

        # Directory to install to
        install_dir = get_python_lib()
        install_dir = self.install_libbase

        # Install files
        print("INSTALL DIR", install_dir)
        # print(self.__dict__)
        [shutil.copy(filename, install_dir) for filename in filenames]

def partition(pred, iterable):
    t1, t2 = itertools.tee(iterable)
    return itertools.filterfalse(pred, t1), filter(pred, t2)


class build(build_orig):

    def finalize_options(self):
        super().finalize_options()
        condition = lambda el: el[0] == 'build_ext'
        rest, sub_build_ext = partition(condition, self.sub_commands)
        self.sub_commands[:] = list(sub_build_ext) + list(rest)

setup(
    name='gridppOI',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='Spatial interpolation using optimal interpolation',

    # The project's main homepage.
    url='https://github.com/metno/gridppOI',

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

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},
    py_modules=['gridpp'],

    cmdclass={'install': CompiledLibInstall},
    # cmdclass={'build': build},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    # entry_points={
    #     'console_scripts': [
    #         'verif=verif:main',
    #     ],
    # },
)
