Installation
====================

**Last updated:** 4/02/2025

LoopStructural is supported and tested on Python 3.9+ and can be installed on Linux, Windows and Mac. 
We recommend installing LoopStructural into clean python environment. Either using anaconda or python virtual environments. 
There are three ways of installing LoopStructural onto your system:

Installing from pip or conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    pip install LoopStructural
    pip install LoopStructural[all] # to include all optional dependencies


.. code-block::

    conda install -c conda-forge -c loop3d loopstructural
    
    
Compiling LoopStructural from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can install the most recent version of LoopStructural by cloning it from GitHub. 




.. code-block::

    git clone https://github.com/Loop3D/LoopStructural.git
    cd LoopStructural
    pip install -e . # -e installs as an editable package so you can modify the source code and see the changes immediately

Dependencies
~~~~~~~~~~~~

Required dependencies:

* numpy
* pandas
* scipy
* scikit-image
* scikit-learn

Optional dependencies:

* matplotlib, 2D/3D visualisation
* pyvista, 3D visualisation
* surfepy, radial basis interpolation
* map2loop, generation of input datasets from regional Australian maps
* geoh5py, export to gocad hdf5 format
* pyevtk, export to vtk format
* dill, serialisation of python objects
* loopsolver, solving of inequalities
* tqdm, progress bar
