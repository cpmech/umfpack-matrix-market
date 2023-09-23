#!/bin/bash

# This script compares the UMFPACK solution with enforced unsymmetric strategy
# and the two pressurized cylinder files

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
./solve_matrix_market pressurized-cylinder-linear-elastic-symmetric-1.mtx 0

echo
echo
echo
./solve_matrix_market pressurized-cylinder-linear-elastic-symmetric-1.mtx 0 1

echo
echo
echo
./solve_matrix_market pressurized-cylinder-linear-elastic-symmetric-2.mtx 0

echo
echo
echo
./solve_matrix_market pressurized-cylinder-linear-elastic-symmetric-2.mtx 0 1
