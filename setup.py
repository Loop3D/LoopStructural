#setup cython code
from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy
setup(
	name="LoopStructural",
	install_requires=[
	'Cython'
	],
        packages=find_packages(),
	# ext_modules=cythonize("LoopStructural/cython/*.pyx",compiler_directives={"language_level": "3"}),
	include_dirs=[numpy.get_include()],
	include_package_data=True,
	)
