SHELL :=/bin/bash
.PHONY: all build notebook notebookbuild
build:
	python3 setup.py install build_ext --inplace;

dependencies:
	sudo apt-get update  && sudo apt-get install python3 python3-venv pybind11-dev mesa-common-dev mesa-utils libgl1-mesa-dev gcc g++; 
	
venv:
ifeq ("","$(wildcard ${LOOP_ENV})")
	python3 -m venv ${LOOP_ENV}
endif
	. ${LOOP_ENV}/bin/activate

all:
	. ${LOOP_ENV}/bin/activate  && pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace; 

notebook:
	. ${LOOP_ENV}/bin/activate &&	jupyter-notebook --no-browser; 

notebookbuild:
	. ${LOOP_ENV}/bin/activate &&	python3 setup.py install build_ext --inplace &&	jupyter-notebook --no-browser; 
compileexamples:
	. ${LOOP_ENV}/bin/activate && . build_notebook.sh;



