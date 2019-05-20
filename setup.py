#setup cython code
from setuptools import setup
from Cython.Build import cythonize
import numpy
setup(
	name="FME",
	install_requires=[
	'Cython'
	],
        packages=['FME'],
	ext_modules=cythonize("FME/dsi_helper.pyx"),
	include_dirs=[numpy.get_include()],
	)
