#setup cython code
from setuptools import setup

setup(
	name = "FME",
        version 
	ext_modules = cythonize("pycompass/SNE/pdf.pyx"),
	include_dirs = [np.get_include(),sp.get_include()]
	)
