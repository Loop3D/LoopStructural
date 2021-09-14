#setup cython code
import sys
try:
	from setuptools import setup, find_packages
except:
	raise RuntimeError("Cannot import setuptools \n"\
	     "python -m pip install setuptools")
	sys.exit(1)

try:
	from Cython.Build import cythonize
except:
	raise RuntimeError("Cannot import cython \n"\
	     "python -m pip install cython")
	sys.exit(1)
try:
	import numpy
except:
	raise RuntimeError("Cannot import numpy \n"\
	     "python -m pip install numpy")
	sys.exit(1)	 
		 
import numpy
import os

package_root = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(package_root, "LoopStructural/version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

setup(
	name="LoopStructural",
	install_requires=[
	'numpy>=1.18', #need to fix numpy to 1.18 because we build against it
	'pandas',
	'scipy',
	'scikit-image',
	'scikit-learn',
	'tqdm',
	],
	version=version,
    packages=find_packages(),
	ext_modules=cythonize("LoopStructural/interpolators/cython/*.pyx",compiler_directives={"language_level": "3"}),
	include_dirs=[numpy.get_include()],
	include_package_data=True,
	package_data={'LoopStructural':['datasets/data/*.csv','datasets/data/*.txt']},
	)
