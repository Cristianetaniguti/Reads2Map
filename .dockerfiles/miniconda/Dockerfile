FROM drtools/alpine-miniconda:4.5.12

RUN conda install -c bioconda -y pyfasta pandas pysam whatshap \
    && conda clean --all
    
COPY ./pedsim_files.py /opt/scripts/

ENTRYPOINT ["/sbin/tini", "--"]