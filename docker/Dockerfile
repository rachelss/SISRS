FROM ubuntu:14.04

ENV SAMTOOLS_VERSION 1.3.1
ENV BOWTIE2_VERSION 2.3.3.1
ENV BBMAP_VERSION 37.66
ENV BAMUTIL_VERSION 1.0.13

# Note: dependencies are grouped below. Since I can't do inline comments to say
# which deps are for what, I'll list them in order here:
# 1. general 2. samtools
RUN apt-get update && apt-get install -y \
    build-essential \
    bzip2 \
    default-jre \
    wget \
    \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    \
    unzip \
  && rm -rf /var/lib/apt/lists/*


RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make install \
    && cd .. \
    && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2
    
# install bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
RUN unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
    && rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
ENV PATH="${PATH}:/bowtie2-${BOWTIE2_VERSION}-linux-x86_64"

# install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz
RUN tar xf BBMap_${BBMAP_VERSION}.tar.gz \
    && rm BBMap_${BBMAP_VERSION}.tar.gz
ENV PATH="${PATH}:${PWD}/bbmap"

# install bamUtil (for bam diff)
RUN wget https://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.${BAMUTIL_VERSION}.tgz
RUN tar xf BamUtilLibStatGen.${BAMUTIL_VERSION}.tgz
RUN cd bamUtil_${BAMUTIL_VERSION} \
    && make all \
    && make install \
    && rm -rf BamUtilLibStatGen.${BAMUTIL_VERSION}.tgz bamUtil_${BAMUTIL_VERSION}

ENV VELVET_VERSION 1.2.10
RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_${VELVET_VERSION}.tgz \
    && tar xf velvet_${VELVET_VERSION}.tgz \
    && cd velvet_${VELVET_VERSION} \
    && make \
    && cp velvet* /usr/local/bin \
    && cd .. \
    && rm -rf velvet_${VELVET_VERSION} velvet_${VELVET_VERSION}.tgz

# install custom python using Miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
RUN bash miniconda.sh -b -p miniconda
ENV PATH="/miniconda/bin:$PATH"
RUN hash -r \
    && conda config --set always_yes yes --set changeps1 no \
    && conda update -q conda \
    && conda info -a \
    && rm miniconda.sh

RUN pip install \
    pytest \
    Biopython

# TODO: move to top for next full rebuild
RUN apt-get update && apt-get install -y \
    parallel \
    libssl-dev \
    git \
  && rm -rf /var/lib/apt/lists/*

# TODO: move up for next full rebuild
RUN git clone https://github.com/DecodeGenetics/BamHash \
    && cd BamHash \
    && make \
    && cp bamhash_checksum_bam /usr/local/bin \
    && cd .. \
    && rm -rf BamHash

COPY run_tests.sh /
