FROM ubuntu:18.04

RUN apt update \
    && apt install -y git build-essential vcflib \
    && git clone --recursive --branch v1.3.0 https://github.com/ekg/freebayes.git \
    && cd freebayes \
    && make