#!/bin/bash

## This script is run by the script from the ../2_analysis/2_compile_facs.ipynb

cd ../1_inputs

curl -O https://zenodo.org/records/15000710/files/raw_FACS.tar.gz

mkdir -p raw_inputs

tar -xvf raw_FACS.tar.gz -C raw_inputs

rm raw_FACS.tar.gz