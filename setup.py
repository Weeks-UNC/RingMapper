#########################################################
# Code to compile readMutStrings.pyx
# Anthony Mustoe
# 2018
#
# This file is licensed under the terms of the MIT license
#
#########################################################
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext = Extension(
    "readMutStrings", sources=["readMutStrings.pyx"], include_dirs=[numpy.get_include()]
)
ext = cythonize(ext, compiler_directives={"language_level": "2"})

setup(
    name="RingMapper",
    ext_modules=ext,
    scripts=["ringmapper.py", "pairmapper.py"],
)
