try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import os

do_cython = True
try:
    from Cython.Build import cythonize
except ImportError:
    do_cython = False

inc = ["../code/",]
try:
    import numpy as np
    inc += [np.get_include(),]
except ImportError:
    do_cython = False

src = [os.path.join("..", "code", fn) for fn in [
    "cyl_point.cpp",
    "point.cpp",
    "GEM.cpp",
    "ICPC.cpp",
    "PPC.cpp",
    "Setup.cpp",
    "Utils.cpp",
    "VelocityLookup.cpp",
    "VelocityModel.cpp"
]]

ext = ".pyx" if do_cython else ".cpp"
src +=  [os.path.join("siggen", "_siggen" + ext)]

extensions = [
    Extension("siggen._siggen", sources=src,
    include_dirs=inc,
    language="c++",
    extra_compile_args=["-std=c++11",
                        "-Wno-unused-function",
                        "-Wno-uninitialized",
                        "-DNO_THREADS"],
    extra_link_args=["-std=c++11"]
    )
]

if do_cython:
    extensions =  cythonize(extensions)

setup(
    name = "siggen",
    version="0.0.1",
    author="Ben Shanks",
    author_email="benjamin.shanks@gmail.com",
    ext_modules =extensions,
    packages=["siggen"],
    install_requires=["numpy", "scipy"]
)
