SHELL :=/bin/bash
.PHONY: all build notebook notebookbuild
build:
	python3 setup.py install build_ext --inplace;

# Ubuntu dependencies
dependencies:
	sudo apt-get update  && sudo apt-get install python3 python3-dev python3-venv pybind11-dev mesa-common-dev mesa-utils libgl1-mesa-dev gcc g++;

# Fedora dependencies
dependencies.fc:
	sudo dnf update  && sudo dnf install python3 python3-devel pybind11-devel mesa-libGL-devel gcc g++;

	
venv:
ifeq ("","$(wildcard ${LOOP_ENV})")
	python3 -m venv ${LOOP_ENV}
endif
	. ${LOOP_ENV}/bin/activate

all:
	pip3 install -r requirements.txt && python3 setup.py install build_ext --inplace;

notebook:
	jupyter-notebook --no-browser;

notebookbuild:
	python3 setup.py install build_ext --inplace &&	jupyter-notebook --no-browser;
compileexamples:
	./build_notebook.sh;



