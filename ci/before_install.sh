#!/bin/bash

# install dependencies available as Ubuntu packages
sudo add-apt-repository --yes ppa:debian-med/ppa
sudo apt-get update -qq
sudo apt-get install -qq samtools bowtie2 parallel mafft

# install BBMap
wget https://sourceforge.net/projects/bbmap/files/BBMap_37.66.tar.gz
tar xzvf BBMap_37.66.tar.gz

# get bamUtil (for diffing bams)
git clone https://github.com/anderspitman/bamUtil_binary

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
