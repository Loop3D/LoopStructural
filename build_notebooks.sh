#!/bin/bash
for filename in ./examples/*.py
	do p2j -t ${filename%.*}_nb.ipynb -o $filename
done
