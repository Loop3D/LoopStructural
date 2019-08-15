SHELL :=/bin/bash
.PHONY: all build notebook notebookbuild

build:
	. ~/fme/bin/activate &&	python3 setup.py install build_ext --inplace;

all:
	sudo apt-get install pybind11-dev mesa-common-dev mesa-utils libl1mesa-dev; 
	. ~/fme/bin/activate  && pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace; 

notebook:
	. ~/fme/bin/activate &&	jupyter-notebook --no-browser; 

notebookbuild:
	. ~/fme/bin/activate &&	python3 setup.py install build_ext --inplace &&	jupyter-notebook --no-browser; 
