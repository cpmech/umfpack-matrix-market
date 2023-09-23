#!/bin/bash

set -e

SOURCE=`pwd`

bash all.bash
cd /tmp/build-umfpack-matrix-market
valgrind --leak-check=full \
    --show-leak-kinds=all \
    --suppressions=${SOURCE}/valgrind.supp \
    ./solve_matrix_market
