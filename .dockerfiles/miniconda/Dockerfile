FROM drtools/alpine-miniconda

RUN conda install -c bioconda -y pyfasta pandas \
    && conda clean --all

COPY ./pedsim_files.py /opt/scripts/