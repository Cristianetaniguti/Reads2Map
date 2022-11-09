FROM ubuntu:18.04

RUN apt update \
    && apt install -y  build-essential wget zlib1g-dev   \
    && wget http://catchenlab.life.illinois.edu/stacks/source/stacks-2.62.tar.gz \
    && tar xfvz stacks-2.62.tar.gz \
    && cd stacks-2.62 \
    && ./configure \
    && make \
    && make install