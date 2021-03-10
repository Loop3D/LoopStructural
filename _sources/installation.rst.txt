Installation
====================
LoopStructural is supported and tested on Python 3.6+ and can be installed on Linux, Windows and Mac. 
We recommend installing LoopStructural into clean python environment. Either using anaconda or python virtual environments. 
There are three ways of installing LoopStructural onto your system:

Installing from pip or conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    pip install LoopStructural

.. code-block::
    conda install -c conda-forge -c loop3d loopstructural

    
Compiling LoopStructural from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can install the most recent version of LoopStructural by cloning it from GitHub. 
You will need to have a C/C++ development environment for compiling cython extensions.

If you are using a linux system you may need to install some dependencies for LavaVu.

.. code-block::

    sudo apt-get update  && sudo apt-get install python3 python3-venv python3-dev make pybind11-dev mesa-common-dev mesa-utils libgl1-mesa-dev gcc g++



.. code-block::

    git clone https://github.com/Loop3D/LoopStructural.git
    cd LoopStructural
    pip install .

Dependencies
~~~~~~~~~~~~

Required dependencies:

* numpy
* pandas
* scipy
* scikit-image
* scikit-learn

Optional dependencies:
* matplotlib, 2D visualisation
* LavaVu, 3D visualisation
* surfepy, radial basis interpolation
* rasterio, exporting triangulated surfaces
* map2loop, generation of input datasets from regional Australian maps


Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LoopStructural can be used either by compiling the docker image or by pulling the compiled
docker image from docker hub.

.. code-block::

    docker pull loop3d/loopstructural
    docker run -i -t -p 8888:8888 -v LOCALDIRPATH:/home/jovyan/shared_volume loop3d/loopstructural`.
