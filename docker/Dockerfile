FROM ubuntu

ARG FASTVIFIVER="1.0"
ARG FASTVIFIUSER="fastvifi"
ARG FASTVIFIGROUP="fastvifi_group"


LABEL dockerfile.version="1.0"
LABEL software="FastViFi"
LABEL software.version="1"
LABEL description="Tool to detect viral presence and hybrid DNA/RNA junctions from NGS data."
LABEL website="https://github.com/sara-javadzadeh/FastViFi"
LABEL license="https://github.com/sara-javadzadeh/FastViFi/blob/main/LICENSE"
LABEL maintainer="Sara Javadzadeh"
LABEL maintainer.email="saraj@eng.ucsd.edu"


RUN groupadd $FASTVIFIGROUP && \
 useradd --home-dir /home -g $FASTVIFIGROUP $FASTVIFIUSER && \
 chown -R ${FASTVIFIUSER}:${FASTVIFIGROUP} /home

USER $FASTVIFIUSER

RUN mkdir /home/fastvifi && \
 mkdir /home/input && \
 mkdir /home/output

USER root

ENV TZ=USA/NY
ENV DEBIAN_FRONTEND=noninteractive
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
#Download libraries for AA
RUN apt-get update
RUN apt-get install -y
RUN apt-get install -y python-dev gfortran python3-pip zlib1g-dev libbz2-dev liblzma-dev libjpeg-dev
RUN apt-get install -y python3.7
RUN update-alternatives --install /usr/local/bin/python python /usr/bin/python3.7 1;
RUN update-alternatives --install /usr/local/bin/python python /usr/bin/python3.6 2;
RUN pip3 install pysam google-cloud
RUN pip3 install numpy scipy matplotlib
RUN apt-get install bwa
RUN apt-get install -y make
RUN apt-get install -y perl
RUN apt-get install -y cpanminus

#Set up gcloud
RUN apt-get install lsb-release curl -y
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)"
RUN echo "deb http://packages.cloud.google.com/apt cloud-sdk-stretch main" >  /etc/apt/sources.list.d/google-cloud-sdk.list
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN apt-get update && \
 apt-get -y install google-cloud-sdk \
 gcc python-dev python-setuptools \
 libncurses5-dev libncursesw5-dev \
 r-base r-base-dev libcurl4-openssl-dev pkg-config \
 git
RUN pip3 install -U crcmod
RUN python3.7 -m pip install --upgrade google-api-python-client
RUN pip3 install requests
RUN pip3 install dendropy

# ViFi uses pysam library on python2.
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python2.7 get-pip.py
RUN python2.7 -m pip install pysam

# Perl Libs
RUN cpanm PerlIO::gzip
RUN apt-get install time
ENV AA_DATA_REPO="/home/data_repo"
ENV PATH="/home/bin:${PATH}"
ENV PYTHONPATH="/home/ViFi/scripts"
ENV VIFI_DIR="/home/ViFi"
ENV REFERENCE_REPO="/home/repo/data"
ENV PERL5LIB=/home/lib/:$PERL5LIB

ADD lib/  /home/lib/
ADD bin/ /home/bin
# Copy compressed versions of some tools from the local machine.
ADD tools/ /home/tools/

WORKDIR /home/
# Install dependencies from source.
RUN tar -jxf /home/tools/bcftools-1.6.tar.bz2 && \
 tar -jxf /home/tools/htslib-1.7.tar.bz2 &&\
 tar -jxf /home/tools/samtools-1.7.tar.bz2
WORKDIR /home/bcftools-1.6
RUN ./configure &&\
 make &&\
 make install
WORKDIR /home/htslib-1.7
RUN ./configure && \
 make &&\
 make install
WORKDIR /home/samtools-1.7
RUN ./configure &&\
 make && \
 make install

WORKDIR /home/

USER $FASTVIFIUSER

# Download ViFi.
RUN cd /home &&\
 git clone https://github.com/sara-javadzadeh/ViFi.git


USER root

# Install dependencies and cleanup apt.
RUN apt-get update && apt-get -y --no-install-recommends install \
 wget \
 ca-certificates \
 zlib1g-dev \
 make \
 g++ \
 rsync \
 git \
 cpanminus && \
 rm -rf /var/lib/apt/lists/* && apt-get autoclean

USER $FASTVIFIUSER
# Perl module required for kraken2-build.
RUN cpanm Getopt::Std

USER $FASTVIFIUSER
# Download Kraken from GitHub and Install.
RUN cd /home && \
 git clone https://github.com/sara-javadzadeh/kraken2.git && \
 cd /home/kraken2 && \
 ./install_kraken2.sh . && \
 mkdir /home/data /home/kraken2-db

ENV PATH="$PATH:/kraken2" \
    LC_ALL=C

# Download FastViFi.
RUN cd /home/fastvifi && \
 git init && \
 git remote add origin https://github.com/sara-javadzadeh/FastViFi.git && \
 git fetch && \
 git checkout main

USER root

RUN chown --changes --no-dereference --recursive \
          ${FASTVIFIUSER}:${FASTVIFIGROUP} \
        /home/fastvifi \
        /home/input \
        /home/data \
        /home/kraken2 \
        /home/kraken2-db


USER $FASTVIFIUSER

WORKDIR /home/fastvifi/FastViFi
