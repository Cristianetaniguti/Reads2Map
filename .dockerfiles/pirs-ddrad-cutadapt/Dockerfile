FROM ubuntu:18.04

# Download pirs and dependencies
RUN apt-get update \
    && apt-get install -y git build-essential libz-dev libboost-all-dev \
    && git clone --branch v2.0.2 https://github.com/galaxy001/pirs.git \
    && cd pirs \
    && make 

# Download ddRADseqTools and dependencies
RUN apt install -y python3-pip \
    && pip3 install numpy \
    && git clone https://github.com/GGFHF/ddRADseqTools.git \
    && cd ddRADseqTools/Package \
    && chmod 777 *.py

# Download cutadapt
RUN apt install -y cutadapt

# Adding pirs library
RUN cd / \
    && mkdir /usr/local/share/pirs/ \
    && mkdir /usr/local/share/pirs/Base-Calling_Profiles  \
    && mkdir /usr/local/share/pirs/InDel_Profiles  \
    && mkdir /usr/local/share/pirs/GC-depth_Profiles  \
    && cp pirs/Profiles/Base-Calling_Profiles/* /usr/local/share/pirs/Base-Calling_Profiles/ \
    && cp pirs/Profiles/InDel_Profiles/* /usr/local/share/pirs/InDel_Profiles/ \
    && cp pirs/Profiles/GC-depth_Profiles/* /usr/local/share/pirs/GC-depth_Profiles/

# Download cleanFastq and dependencies
RUN git clone https://github.com/davidvi/cleanFastq.git \
    && cd cleanFastq \
    && g++ -std=c++11 fixFastq.cpp -o fixFastq


WORKDIR /