#!/bin/bash

# Params to change are P Q in seq P Q
# Where P--Q is the range of county rasters and X is the number of cores to use

seq 1 410 | parallel --no-notice -j 8 python ./SOC_LEM_County.py
