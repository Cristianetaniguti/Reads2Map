FROM cristaniguti/r-samtools:latest

COPY ./simuscopR_0.0.0.9000.tar.gz /opt/

RUN R -e "install.packages('R.utils', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('devtools', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('vcfR', repos='http://brieger.esalq.usp.br/CRAN/');\
          install.packages('tidyverse', repos='http://brieger.esalq.usp.br/CRAN/')"
          
RUN R -e "install.packages('/opt/simuscopR_0.0.0.9000.tar.gz', type = 'source', repos=NULL)"