FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
    git=1:2.25.1-1ubuntu3 \
    python3-pip \
    tabix \
    libbz2-dev \
    liblzma-dev

RUN python3 -m pip install numpy==1.22 scipy numba pysam pytest # requirements.txt

RUN git clone --branch v0.6.0 https://github.com/PlantandFoodResearch/MCHap.git   \
    && cd MCHap \
    && python3 setup.py sdist \
    && python3 -m pip install dist/mchap-*tar.gz
