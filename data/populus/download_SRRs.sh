for i in $(cat SRRs.txt); do

    echo $i
    docker run -v $(pwd):/opt/ cristaniguti/sratoolkit ./fasterq-dump $i -O /opt/

done
