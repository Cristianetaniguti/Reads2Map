FROM drtools/alpine-miniconda

RUN conda install -c bioconda -y pyfasta pandas pysam \
    && conda clean --all

COPY ./pedsim_files.py /opt/scripts/