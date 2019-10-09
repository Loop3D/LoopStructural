# LoopStructural
Loop3D Geological Forward Modeling Engine.

* Python/cython implementation of a piecewise linear interpolation (Discrete Smooth Interpolator in Gocad) 
* Fold interpolation using constraints outlined in Laurent 2016 with fold geostatistical tools shown in Grose et al., 2017
* Fault interpolation 

If you want to use LoopStructural the best way to use it is within a jupyter notebook environment. 

## Using docker
LoopStructural can then be used by first building the docker image `docker build -t=LoopStructural .`.
LoopStructural can be used by run  `docker run  -i -t -p 8888:8888 -v $WORKDIR:/notebooks LoopStructural ` where **`$WORKDIR`** is the path to your jupyter notebooks. You must change this to the absolute path! This will start a notebook server running on localhost:8888 without  password or certificate. Be aware any changes made to the notebooks within the notebooks directory will be saved on your local versions.

## On Linux
LoopStructural can be easily installed using a Makefile once a few things are set up. Firstly, you need to add an environment variable to your system. LoopStructural_ENV, this can be done by adding `export LOOP_VENV=$YOUR_PATH_TO_VIRTUAL_ENVIRONMENT` to the `.bashrc` file. Make sure the path is updated to a directory in your system where you want to save the python virtual environment. It could be for example where you clone this repository and a subfolder called venv or LoopStructural. 

Once you have the environment variable you can run the command `make dependencies` which will install the required dependencies for LoopStructural:

Required dependencies
* pybind11-dev
* mesa-common-dev
* mesa-utils
* libgl1mesa-dev
`sudo apt-get install pybind11-dev mesa-common-dev mesa-utils libl1mesa-dev`

Once these are installed you can run `make venv` to create a new python virtual environment in the location you specified. If a python environment already exists then this will be used.

`make all` will install the required python dependencies for LoopStructural and then install and build the library. It just executes the following command: `pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace`

If you want to use a jupyter notebook then you can launch a server by running `make notebook`, alternatively you can run `make notebookbuild` if you want to build the library before launching the server.

If you want to compile the example files into jupyter notebooks you can do this using the `p2j` package. This can be done by running `make compileexamples`

## On windows
You can install LoopStructural on windows natively if you have a working C compiler. Otherwise the easiest way is to use the linux subsystem for linux or docker. 
In order to use lavavu you will have to use some method of fowarding X to the windows client. 

## Problems
Any bugs/feature requests/comments send to lachlan.grose@monash.edu
