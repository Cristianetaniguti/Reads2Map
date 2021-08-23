FROM ubuntu:20.04

RUN apt-get update \
    && apt-get install -y \
    gfortran=4:9.3.0-1ubuntu2 \
    libopenblas-dev=0.3.8+ds-1ubuntu0.20.04.1 \
    liblapack-dev=3.9.0-1build1 \
    libgsl-dev=2.5+dfsg-6build1 \
    git=1:2.25.1-1ubuntu3 \
    vcftools=0.1.16-1build1 

RUN apt-get install -y \
    python3-pip \
    && python3 -m pip install radinitio 

RUN git clone https://github.com/lh3/seqtk.git \
    && cd seqtk \
    && make

ENV PATH="/seqtk:${PATH}"       

CMD bash