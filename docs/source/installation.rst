FME Installation
================
FME is a python3 library with a few dependencies. A makefile exists to make the installation process easier on unix systems. 
For a windows use it is recommended to use the `Windows Linux Subsystem <https://docs.microsoft.com/en-us/windows/wsl/install-win10>` for running WLS. 
**Native installation on windows is not recommended** because it can be complicated due to windows missing a C++ compiler. Final releases of FME will provide compiled wheels that can be installed. 

Dependencies
------------
* `meshpy <https://documen.tician.de/meshpy/installation.html>` is used to build tetrahedral meshes
* `lavavu <https://github.com/lavavu/LavaVu>` is used for 3D visualisation
* numpy
* scipy
* scikit-learn
* scikit-image
* meshio
* cython
  
Other useful packages that are used in our examples are:
* matplotlib
* geopandas 
Docker
------
A docker file is 

Windows
--------

Linux
-----

FME can be easily installed using a Makefile once a few things are set up. Firstly, you need to add an environment variable to your system. FME_ENV, this can be done by adding `export FME_VENV=$YOUR_PATH_TO_VIRTUAL_ENVIRONMENT` to the `.bashrc` file. Make sure the path is updated to a directory in your system where you want to save the python virtual environment. It could be for example where you clone this repository and a subfolder called venv or fme. 

Once you have the environment variable you can run the command `make dependencies` which will install the required dependencies for FME:

Required dependencies
* pybind11-dev
* mesa-common-dev
* mesa-utils
* libgl1mesa-dev
`sudo apt-get install pybind11-dev mesa-common-dev mesa-utils libl1mesa-dev`

Once these are installed you can run `make venv` to create a new python virtual environment in the location you specified. If a python environment already exists then this will be used.

`make all` will install the required python dependencies for FME and then install and build the library. It just executes the following command: `pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace`

If you want to use a jupyter notebook then you can launch a server by running `make notebook`, alternatively you can run `make notebookbuild` if you want to build the library before launching the server.

If you want to compile the example files into jupyter notebooks you can do this using the `p2j` package. This can be done by running `make compileexamples`

