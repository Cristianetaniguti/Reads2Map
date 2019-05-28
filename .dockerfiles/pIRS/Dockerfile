FROM ubuntu:18.04

# Download pirs and dependencies
RUN apt-get update \
    && apt-get install -y git build-essential libz-dev libboost-all-dev \
    && git clone --branch v2.0.2 https://github.com/galaxy001/pirs.git \
    && cd pirs \
    && make

WORKDIR pirs/src/pirs/
