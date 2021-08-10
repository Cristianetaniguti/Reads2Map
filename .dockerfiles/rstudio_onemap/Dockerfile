FROM cristaniguti/onemap_git:latest

# For supermassa
RUN apt-get update \
    && apt install -y python2 python-tk curl

RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py 
RUN python2 get-pip.py 
RUN pip install matplotlib==2.0.2 numpy 
RUN ln -s /usr/bin/python2 /usr/bin/python

RUN R -e "install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('stringi', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('doParallel', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('matrixStats', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('BiocManager', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shiny', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shinydashboard', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('shinymanager', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('tidyverse', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('largeList', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('ggpubr', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('config', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('irr', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('RColorBrewer', repos='http://brieger.esalq.usp.br/CRAN/'); \
          install.packages('polyRAD', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('pROC', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('DT', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('golem', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('htmltools', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('htmlwidgets', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vroom', repos='http://brieger.esalq.usp.br/CRAN/');\
          BiocManager::install('pcaMethods');\
          BiocManager::install('VariantAnnotation');\
          devtools::install_github('tpbilton/GUSbase');\
          devtools::install_github('tpbilton/GUSMap')"

COPY ./onemapUTILS_0.2.0.tar.gz /opt/
RUN R -e "install.packages('/opt/onemapUTILS_0.2.0.tar.gz', type = 'source', repos=NULL)"

RUN Rscript -e "install.packages('filelock')"
RUN Rscript -e "install.packages('vroom')"
