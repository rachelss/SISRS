#!/bin/bash

# install dependencies available as Ubuntu packages
sudo add-apt-repository --yes ppa:debian-med/ppa
sudo apt-get update -qq
sudo apt-get install -qq samtools bowtie2 python-biopython parallel mafft

# install BBMap
wget https://sourceforge.net/projects/bbmap/files/BBMap_37.66.tar.gz
tar xzvf BBMap_37.66.tar.gz
export PATH=$PATH:~/bbmap/

# download test data
git clone git@github.com:BobLiterman/SISRS_Small_Data
