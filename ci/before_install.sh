#!/bin/bash

sudo add-apt-repository --yes ppa:debian-med/ppa
sudo apt-get update -qq
sudo apt-get install -qq samtools bowtie2 python-biopython parallel mafft
wget https://sourceforge.net/projects/bbmap/files/BBMap_37.56.tar.gz
tar xzvf BBMap_37.56.tar.gz
export PATH=$PATH:~/bbmap/


