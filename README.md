# fast-kmedoids

*Copyright (c) 2021-2022 by Miguel A. Caro and Patricia Hernández-León*

**fast-kmedoids** is a kmedoids Python implementation that uses Fortran libraries and
f2py for faster execution. It has been written by **Miguel A. Caro** and **Patricia
Hernández-León** based on Cristian Bauckhage's implementation.

**fast-kmedoids** is released under the GNU General Public License version 3, see
LICENSE.md file.

## Installation

### Prerrequisites

- Numpy
- A Fortran compiler (successfully tested with `gfortran`)

### Building the libraries

Execute the build script:

    ./build_libraries.sh

Add the root directory to your Python path:

    echo "export PYTHONPATH=$(pwd):\$PYTHONPATH" >> ~/.bashrc
    source ~/.bashrc

## Usage

See `example.py` for a typical use case.
