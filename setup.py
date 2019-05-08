#setup cython code
from setuptools import setup
from Cython.Build import cythonize

setup(
	name="FME",
        packages=['FME'],
	ext_modules=cythonize("FME/dsi_helper.pyx"),
	)
