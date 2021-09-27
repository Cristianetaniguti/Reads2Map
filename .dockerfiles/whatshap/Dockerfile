FROM conda/miniconda3:latest

RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge

RUN conda install whatshap
