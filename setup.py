#########################################################
# Code to compile readMutStrings.pyx
# Anthony Mustoe
# 2018
#
# This file is licensed under the terms of the MIT license
#
#########################################################
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

ext = Extension(
    name="ringmapper.readMutStrings",
    sources=["ringmapper/readMutStrings.pyx"],
    include_dirs=[numpy.get_include()],
)
ext = cythonize(ext, compiler_directives={"language_level": "2"})

setup(
    name="RingMapper",
    packages=["ringmapper"],
    package_dir={"ringmapper": "./"},
    package_data={"ringmapper": ["*.pxd"]},
    # include_package_data=True,
    py_modules=[
        "ringmapper.ringmapper",
        "ringmapper.pairmapper",
        "ringmapper.ShapeMapper_MMS0_Mut_Filter",
        "ringmapper.readMutStrings",
    ],
    ext_modules=ext,
    scripts=[
        "./ringmapper/ringmapper.py",
        "./ringmapper/pairmapper.py",
        "./ringmapper/ShapeMapper_MMS0_Mut_Filter.py",
    ],
)
