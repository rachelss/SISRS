#!/bin/bash

SAMTOOLS_VERSION=1.3.1
BOWTIE2_VERSION=2.3.3.1
BBMAP_VERSION=37.66
SISRS_DIR=$PWD

# install dependencies available as Ubuntu packages
apt-get update -y
apt-get install -y wget build-essential bzip2 \
    parallel mafft default-jre

touch $SISRS_DIR/ci/set_env.sh

# build and install samtools
cd $HOME
apt-get install -y liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev 
wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}
./configure
make install
cd ..
cd $SISRS_DIR

# install bowtie2
cd $HOME
wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
echo export PATH='$PATH':$PWD/bowtie2-${BOWTIE2_VERSION}-linux-x86_64 >> $SISRS_DIR/ci/set_env.sh
cd $SISRS_DIR

# install BBMap
wget https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz
tar xzvf BBMap_${BBMAP_VERSION}.tar.gz

# get bamUtil (for diffing bams)
git clone --depth 1 https://github.com/anderspitman/bamUtil_binary

# download test data
#git clone https://github.com/BobLiterman/SISRS_Small_Data
git clone --depth 1 https://github.com/anderspitman/SISRS_ci_data pipeline_stages

# install custom python using Miniconda
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
# Useful for debugging any issues with conda
conda info -a

conda create -q -n sisrs-python-env python=2.7
source activate sisrs-python-env

pip install Biopython
pip install pytest
pip install -e .
