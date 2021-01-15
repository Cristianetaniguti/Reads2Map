FROM ubuntu:20.04

RUN apt-get update \
    && apt-get install -y gfortran libopenblas-dev liblapack-dev libgsl-dev git vcftools 

RUN apt-get install -y python3-pip \
    && python3 -m pip install radinitio 

RUN git clone https://github.com/lh3/seqtk.git \
    && cd seqtk \
    && make

ENV PATH="/seqtk:${PATH}"       

CMD bash