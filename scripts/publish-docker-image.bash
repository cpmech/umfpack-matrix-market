#!/bin/bash

NAME="cpmech/umfpack-matrix-market"
VERSION="latest"

docker logout
docker login
docker push $NAME:$VERSION
