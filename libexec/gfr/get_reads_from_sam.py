#!/usr/bin/env python2
#takes bam file aligned to ref_genes for one species
#outputs .fa file for each locus for one species

import sys
from decimal import *
import pysam
from collections import Counter,defaultdict

infile = sys.argv[1]
reads=defaultdict(list)
mainfolder_name='/'.join(infile.split('/')[:-1])      #everything but sam file name
samfile = open(infile,'r')
for line in samfile:
    splitline = line.split()
    if int(splitline[4])>10:
        reads[splitline[2]].append((splitline[0],splitline[9]))    #note: read is already rev-comp if necessary
    if splitline[6] !='-' and int(splitline[8])>10:     #make sure working with current read secondary mapping not next read
        reads[splitline[2]].append((splitline[0],splitline[9]))       #include secondary mapping if good
samfile.close()

for contig,seqs in reads.iteritems():
    outfile = open(contig+'.fa','w')
    for i in seqs:
        outfile.write('>'+i[0]+"\n"+i[1]+"\n")
outfile.close()