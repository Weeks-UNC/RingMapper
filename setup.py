#########################################################
# Code to compile readMutStrings.pyx
# Anthony Mustoe
# 2018
#
# This file is licensed under the terms of the MIT license
#
#########################################################
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy

ext = Extension(
    name="smccp.readMutStrings",
    sources=["smccp/readMutStrings.pyx"],
    include_dirs=[numpy.get_include()],
)
ext = cythonize(ext, compiler_directives={"language_level": "2"})

setup(
    name="RingMapper",
    packages=find_packages(include=["smccp", "smccp.*"]),
    package_data={"": ["*.pxd"]},
    # include_package_data=True,
    ext_modules=ext,
    scripts=[
        "./smccp/ringmapper.py",
        "./smccp/pairmapper.py",
        "./smccp/ShapeMapper_MMS0_Mut_Filter.py",
    ],
)
