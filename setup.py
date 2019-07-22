#setup cython code
from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy
setup(
	name="FME",
	install_requires=[
	'Cython'
	],
        packages=find_packages(),
	ext_modules=cythonize("FME/cython/*.pyx"),
	include_dirs=[numpy.get_include()],
	)
