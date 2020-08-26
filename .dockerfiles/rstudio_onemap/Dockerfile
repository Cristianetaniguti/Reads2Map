FROM rocker/rstudio:3.6.2

RUN apt update \
   && apt-get install -y libcurl4-gnutls-dev libssl-dev libftgl2 libglu1-mesa-dev libftgl2 libfreetype6-dev  libgfortran5 python-pip python-tk libxml2-dev libjpeg-dev libbz2-dev liblzma-dev qpdf \
   && pip install matplotlib==2.0.2 numpy 

RUN R -e "install.packages('updog', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringi', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('reshape2', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('ggplot2', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('doParallel', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('matrixStats', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('BiocManager', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shiny', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shinydashboard', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('tidyverse', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('largeList', repos='http://brieger.esalq.usp.br/CRAN/'); \
          BiocManager::install('pcaMethods');\
          BiocManager::install('VariantAnnotation' );\
          install.packages('polyRAD', repos='http://brieger.esalq.usp.br/CRAN/');\
          devtools::install_github('tpbilton/GUSbase');\
          devtools::install_github('tpbilton/GUSMap')"

RUN R -e "devtools::install_github('Cristianetaniguti/onemap')"

RUN R -e "devtools::install_github('Cristianetaniguti/genotyping4onemap', auth_token ='1bcc1d0962f96c9b4d9afc936b807ab2d0804d01')"

RUN wget https://cran.r-project.org/src/contrib/Archive/vcfR/vcfR_1.10.0.tar.gz 

RUN apt update \
    && apt install -y devscripts

RUN R -e "install.packages('vcfR_1.10.0.tar.gz', repos=NULL, type='source')"

RUN R -e  "install.packages('ggpubr', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('config', repos='http://brieger.esalq.usp.br/CRAN/')"
