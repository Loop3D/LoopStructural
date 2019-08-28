#!/bin/bash
#. ~/fme/bin/activate;
#pip3 install p2j
. $FME_ENV
for filename in ./examples/*.py;
	do p2j -t ${filename%.*}_nb.ipynb -o $filename;
done
