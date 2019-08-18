#!/bin/bash
#. ~/fme/bin/activate;
#pip3 install p2j
for filename in ./examples/*.py;
	do p2j -t ${filename%.*}_nb.ipynb $filename;
done
