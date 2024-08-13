import sys

import numpy as np
import setuptools
from Cython.Build import cythonize

if sys.platform == "win32":
    openmp_args = ["/openmp"]
else:
    openmp_args = ["-fopenmp"]


ext_modules = [
    setuptools.Extension(
        "*",
        ["**/*.pyx"],
        extra_compile_args=openmp_args,
        extra_link_args=openmp_args,
        include_dirs=[np.get_include()],
    ),
]


setuptools.setup(
    ext_modules=cythonize(ext_modules),
)
