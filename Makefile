SHELL :=/bin/bash
.PHONY: all build notebook notebookbuild
build:
	. ${FME_ENV} &&	python3 setup.py install build_ext --inplace;

dependencies:
	sudo apt-get install python3 python3-venv pybind11-dev mesa-common-dev mesa-utils libgl1-mesa-dev; 
	
venv:
ifeq ("","$(wildcard ${FME_ENV})")
	python3 -m venv ${FME_ENV}
endif
	. ${FME_ENV}

all:
	. ${FME_ENV}  && pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace; 

notebook:
	. ${FME_ENV} &&	jupyter-notebook --no-browser; 

notebookbuild:
	. ${FME_ENV} &&	python3 setup.py install build_ext --inplace &&	jupyter-notebook --no-browser; 
compileexamples:
	. ${FME_ENV} && . build_notebook.sh;



