#!/bin/bash

DIR=$PWD

source ci/set_env.sh
export PATH="$HOME/miniconda/bin:$PATH"
export PATH=$PATH:$DIR/bbmap/
export PATH=$PATH:$DIR/bamUtil_binary

source activate sisrs-python-env
pytest -s

#bin/sisrs outputAlignment -a premade -c 0 -f $DIR/SISRS_Small_Data/2_identifyFixedSites/ -z $DIR/output/
#cmp $DIR/output/alignment.nex $DIR/SISRS_Small_Data/3_outputAlignment/alignment.nex
#cmp $DIR/output/alignment_bi.nex $DIR/SISRS_Small_Data/3_outputAlignment/alignment_bi.nex
#cmp $DIR/output/alignment_pi.nex $DIR/SISRS_Small_Data/3_outputAlignment/alignment_pi.nex
