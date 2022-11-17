#!/usr/bin/env bash
set -e

for dockerfile in */Dockerfile; do
    DIR=$(dirname "$dockerfile")
    cd "$DIR" && docker build -t "$DIR":v1 . && cd ..
done
