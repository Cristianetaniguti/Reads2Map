FROM ubuntu:20.04

RUN apt-get update \
    && apt-get install -y \
    git=1:2.25.1-1ubuntu3 \
    python3-pip \
    tabix \
    && python3 -m pip install numba scipy pysam networkx graphviz pytest # requirements.txt

RUN git clone https://github.com/PlantandFoodResearch/MCHap.git \
    && cd MCHap \
    && python3 setup.py sdist \
    && python3 -m pip install dist/mchap-*tar.gz
