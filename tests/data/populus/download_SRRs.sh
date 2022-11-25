for i in $(cat SRRs.txt); do

    echo $i
    echo docker run -v $(pwd):/opt/ --rm --entrypoint /usr/bin/fastq-dump -w /opt/ cyverseuk/fastq-dump:latest --gzip $i

done
