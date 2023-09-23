#!/bin/bash

set -e

SOURCE=`pwd`

cd /tmp
mkdir -p build-umfpack-matrix-market
cd build-umfpack-matrix-market/
cmake -S $SOURCE 
make && ./solve_matrix_market
