FROM ubuntu:18.04

RUN apt update

# Dependencies
RUN apt install -y git build-essential libz-dev  libboost-all-dev

# Download pirs
RUN git clone https://github.com/galaxy001/pirs.git

# Install pirs
RUN cd pirs && make

WORKDIR pirs/src/pirs/


