#!/usr/bin/env python3
"""
Setup script for QuasiLocal Toolkit Python bindings.
"""

import os
import sys
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Project metadata
NAME = "quasilocal-toolkit"
VERSION = "0.1.0"
DESCRIPTION = "QuasiLocal Toolkit for numerical relativity"
LONG_DESCRIPTION = """
Mathematically accurate implementation of quasilocal conserved quantities
from "Geometric Vorticity and Quasilocal Angular Momentum Conservation"
by Rayan D. Peru.

Provides tools for computing angular momentum, mass, and other conserved
quantities in numerical relativity without horizon finding.
"""
AUTHOR = "Rayan D. Peru"
LICENSE = "MIT"
URL = "https://github.com/username/quasilocal-toolkit"

# CMake extension builder
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))
        
        for ext in self.extensions:
            self.build_extension(ext)
    
    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        # Required for CMake
        extdir = os.path.join(extdir, NAME)
        
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DBUILD_PYTHON_BINDINGS=ON',
            '-DBUILD_TESTS=OFF',  # Don't build tests for pip package
            '-DBUILD_EXAMPLES=OFF',
        ]
        
        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ['--', '-j2']
        
        env = os.environ.copy()
        
        # Build directory
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        
        # Configure
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                             cwd=self.build_temp, env=env)
        
        # Build
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                             cwd=self.build_temp)

# Setup
setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author=AUTHOR,
    license=LICENSE,
    url=URL,
    
    # Python package structure
    packages=['qlt'],
    package_dir={'qlt': '.'},
    
    # CMake extension
    ext_modules=[CMakeExtension('qlt')],
    cmdclass={'build_ext': CMakeBuild},
    
    # Dependencies
    install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.0',
        'h5py>=3.0.0',
        'matplotlib>=3.3.0',
    ],
    
    # Optional dependencies
    extras_require={
        'docs': ['sphinx>=4.0', 'numpydoc'],
        'tests': ['pytest>=6.0', 'pytest-benchmark'],
        'full': ['mpi4py', 'pandas', 'xarray'],
    },
    
    # Classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    
    # Keywords
    keywords='numerical-relativity general-relativity black-holes angular-momentum',
    
    # Entry points
    entry_points={
        'console_scripts': [
            'qlt-cli=qlt.cli:main',
        ],
    },
    
    # Include package data
    include_package_data=True,
    zip_safe=False,
    
    # Python requirement
    python_requires='>=3.8',
)
