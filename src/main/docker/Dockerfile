FROM ubuntu:14.04
LABEL maintainer "sagory-zalkind@hotmail.fr"


# install packages required for installation
RUN apt-get update && apt-get install -y \
    build-essential \
    python3.4 \
    wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Install NCBI blast+
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz \
    && tar zxf ncbi-blast-2.6.0+-x64-linux.tar.gz \
    && mv ncbi-blast-2.6.0+ /usr/local/ \
    && cd /usr/local/ncbi-blast-2.6.0+/ \
    && mkdir db


# copy python source code
COPY src/main/python/* /


# copy entrypoint
COPY src/main/docker/docker-entrypoint.sh /


# make python code executable
RUN chmod 755 identification16S.py \
    && chmod 755 refineRDPdatabase.py \
    && chmod 755 docker-entrypoint.sh


# copy database to expected folder for blast
ADD src/main/resources/16S_named_RDP.tar.gz /usr/local/ncbi-blast-2.6.0+/db/


# add a volume to read FASTA file
VOLUME /data


# entrypoint (change user specified parameters, initialize database and run pipeline)
ENTRYPOINT ["/docker-entrypoint.sh"]