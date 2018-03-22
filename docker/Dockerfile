FROM ubuntu:14.04

# Note: dependencies are grouped below. Since I can't do inline comments to say
# which deps are for what, I'll list them in order here:
# 1. general 2. samtools 3. BamHash 4. sisrs
RUN apt-get update && apt-get install -y \
    build-essential \
    bzip2 \
    default-jre \
    wget \
    unzip \
    \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    \
    git \
    libssl-dev \
    \
    parallel \
  && rm -rf /var/lib/apt/lists/*


#install samtools
ENV SAMTOOLS_VERSION 1.3.1
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make install \
    && cd .. \
    && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2
    
# install bowtie2
ENV BOWTIE2_VERSION 2.3.3.1
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
RUN unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
    && rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
ENV PATH="${PATH}:/bowtie2-${BOWTIE2_VERSION}-linux-x86_64"

# install BBMap
ENV BBMAP_VERSION 37.66
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz
RUN tar xf BBMap_${BBMAP_VERSION}.tar.gz \
    && rm BBMap_${BBMAP_VERSION}.tar.gz
ENV PATH="${PATH}:${PWD}/bbmap"

# install velvet
ENV VELVET_VERSION 1.2.10
RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_${VELVET_VERSION}.tgz \
    && tar xf velvet_${VELVET_VERSION}.tgz \
    && cd velvet_${VELVET_VERSION} \
    && make \
    && cp velvet* /usr/local/bin \
    && cd .. \
    && rm -rf velvet_${VELVET_VERSION} velvet_${VELVET_VERSION}.tgz

# install mafft
ENV MAFFT_VERSION 7.313
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz
RUN tar xf mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
    && cd mafft-${MAFFT_VERSION}-without-extensions/core \
    && make install \
    && cd ../../ \
    && rm -rf mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
    && rm -rf mafft-${MAFFT_VERSION}-without-extensions

# install BamHash for comparing bam files
RUN git clone https://github.com/DecodeGenetics/BamHash \
    && cd BamHash \
    && git checkout 1c987784929f5c6cf28b9d36d00017380b922c71 \
    && make \
    && cp bamhash_checksum_bam /usr/local/bin \
    && cd .. \
    && rm -rf BamHash

# install custom python using Miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
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

COPY run_tests.sh /
COPY bamdiff /usr/local/bin/
