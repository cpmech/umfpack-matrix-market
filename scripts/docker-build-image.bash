#!/bin/bash

# image name
NAME="cpmech/umfpack-matrix-market"

# build Docker image
docker build --no-cache -t $NAME .

echo
echo
echo
echo "... SUCCESS ..."
echo
echo "... image ${NAME} created ..."
echo
