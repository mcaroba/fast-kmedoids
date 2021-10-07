#!/bin/bash

cd src
python3 -m numpy.f2py --f90flags='-fopenmp' -lgomp -c kmedoids.f90 -m fast_module
