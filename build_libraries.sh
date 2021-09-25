#!/bin/bash

cd src
f2py3 --f90flags='-fopenmp' -lgomp -c kmedoids.f90 -m fast_module
