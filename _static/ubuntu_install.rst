
Linux
~~~~~~~~~~~~
.. container:: toggle

    .. container:: header

        **Development environment**


    LoopStructural can be easily installed using a Makefile once a few things are set up. Firstly, you need to add an
    environment variable to your system. :code:`LOOP_ENV`, this can be done by adding

    .. code-block::

        export LOOP_ENV=$YOUR_PATH_TO_VIRTUAL_ENVIRONMENT

    to your :code:`.bashrc` file.
    Make sure the path is updated to a directory in your system where you want to save the python virtual environment.
    It could be for example where you clone this repository and a subfolder called venv or LoopStructural.

    Once you have the environment variable you can run the command :code:`make dependencies` (or :code:`make dependencies.fc` for Fedora) which will install the required dependencies for LoopStructural:

    Required dependencies for Ubuntu

    * python3
    * python3-dev
    * python3-venv
    * pybind11-dev
    * mesa-common-dev
    * mesa-utils
    * libgl1-mesa-dev
    * gcc
    * g++

    .. code-block::

        sudo apt-get install python3 python3-dev python3-venv pybind11-dev mesa-common-dev mesa-utils libgl1-mesa-dev gcc g++

    Required dependencies for Fedora

    * python3
    * python3-devel
    * pybind11-devel
    * mesa-libGL-devel
    * gcc
    * g++

    .. code-block::

        sudo dnf install python3 python3-devel pybind11-devel mesa-libGL-devel gcc g++

    Once these are installed you can run :code:`make venv` to create a new python virtual environment in the location you
    specified. If a python environment already exists then this will be used.

    :code:`make all` will install the required python dependencies for LoopStructural and then install and build the library.
    It just executes the following command:

    .. code-block::

        pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace

    If you want to use a jupyter notebook then you can launch a server by running :code:`make notebook`, alternatively you can
    run :code:`make notebookbuild` if you want to build the library before launching the server.


