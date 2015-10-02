#!/usr/bin/env python2
#takes bam file aligned to ref_genes for one species
#outputs .fa file for all loci for one species

import sys
from decimal import *
import pysam
from collections import Counter

#######################
def get_consensus(data): 
    pos_bases_counts = Counter(data).most_common()
    if len(pos_bases_counts)>1:
        prop = Decimal(pos_bases_counts[0][1]) / (Decimal(pos_bases_counts[1][1])+Decimal(pos_bases_counts[0][1]))
        if prop > 0.8:
            c_base = pos_bases_counts[0][0]
        else:
            c_base = 'N'
    elif len(pos_bases_counts)==1:
        c_base = pos_bases_counts[0][0]
    else:
        c_base = 'N'
    
    return c_base

def get_consensus_MR(data): 
    pos_bases_counts = Counter(data).most_common()
    if len(pos_bases_counts)>0:
        c_base = pos_bases_counts[0][0]
    else:
        c_base = 'N'
    
    return c_base

######################
infile = sys.argv[1]
minread = int(sys.argv[2])
mainfolder_name='/'.join(infile.split('/')[:-2])      #everything but bam file name and sample name
sp_name=infile.split('/')[-2]      #just sample name
sp_name=sp_name.replace('_loci','')      #just sample name
final_seqs={}
outfile = open(mainfolder_name+'/loci/'+sp_name+'.fa','w')

samfile = pysam.AlignmentFile(infile, "rb" )       #open bamfile
for pileupcolumn in samfile.pileup():       #go through each position for each contig
    contig = samfile.getrname(pileupcolumn.reference_id)
    if contig not in final_seqs:
        final_seqs[contig]=[]
    if int(pileupcolumn.n) >= minread:     #min 3 reads aligned
        bases_aligned = [pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in pileupcolumn.pileups if type(pileupread.query_position) is int]     #this is the pileup
        if len(bases_aligned) >= minread:
            consensus_base = get_consensus_MR(bases_aligned)
        else:
            consensus_base = 'N'
    else:
        consensus_base = 'N'
    final_seqs[contig].append(consensus_base)

for contig,seq in final_seqs.iteritems():
    outfile.write('>'+contig+"\n")
    outfile.write(''.join(seq)+"\n")

outfile.close()