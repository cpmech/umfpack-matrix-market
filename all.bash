#!/bin/bash

set -e

SOURCE=`pwd`

cd /tmp
mkdir -p build-umfpack-matrix-market
cd build-umfpack-matrix-market/
cmake -S $SOURCE 
make

echo
echo
echo
./solve_matrix_market bfwb62.mtx 1 0

echo
echo
echo
./solve_matrix_market bfwb62.mtx 1 1
