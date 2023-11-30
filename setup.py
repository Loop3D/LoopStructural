# setup cython code
import sys

try:
    from setuptools import setup, find_packages
except:
    raise RuntimeError("Cannot import setuptools \n" "python -m pip install setuptools")
    sys.exit(1)

try:
    from Cython.Build import cythonize
except:
    raise RuntimeError("Cannot import cython \n" "python -m pip install cython")
    sys.exit(1)
try:
    import numpy
except:
    raise RuntimeError("Cannot import numpy \n" "python -m pip install numpy")
    sys.exit(1)

import numpy
import os
import codecs

package_root = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(package_root, "LoopStructural/version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

setup(
    name="LoopStructural",
    install_requires=[
        "numpy>=1.18",
        "pandas",
        "scipy",
        "scikit-image",
        "scikit-learn",
        "tqdm",
    ],
    description="Open source 3D structural geology modelling",
    long_description=codecs.open("README.md", "r", "utf-8").read(),
    author="Lachlan Grose",
    author_email="lachlan.grose@monash.edu",
    license=("MIT"),
    url="https://loop3d.github.io/LoopStructural/",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Other Audience",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "License :: Free for non-commercial use",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Multimedia :: Graphics :: 3D Modeling",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: GIS",
    ],
    version=version,
    packages=find_packages(),
    ext_modules=cythonize(
        "LoopStructural/interpolators/_cython/*.pyx",
        compiler_directives={"language_level": "3"},
    ),
    include_dirs=[numpy.get_include()],
    include_package_data=True,
    package_data={
        "LoopStructural": [
            "datasets/data/fault_trace/*",
            "datasets/data/*.csv",
            "datasets/data/*.txt",
            "datasets/data/geological_map_data/*.csv",
            "datasets/data/geological_map_data/*.txt",
        ]
    },
    keywords=[
        "earth sciences",
        "geology",
        "3-D modelling",
        "structural geology",
        "uncertainty",
    ],
)
