#!/bin/bash

export PATH=$PATH:~/bbmap/

bin/sisrs outputAlignment -a premade -c 0 -f ~/SISRS_Small_Data/2_identifyFixedSites/ -z output/
cmp output/3_outputAlignment/alignment.nex SISRS_Small_Data/3_outputAlignment/alignment.nex
cmp output/3_outputAlignment/alignmentbi.nex SISRS_Small_Data/3_outputAlignment/alignmentbi.nex
cmp output/3_outputAlignment/alignmentpi.nex SISRS_Small_Data/3_outputAlignment/alignmentpi.nex
