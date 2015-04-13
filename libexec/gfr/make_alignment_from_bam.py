#!/usr/bin/env python2
#takes bam file aligned to ref_genes
#outputs .fa file for each locus

import sys
from decimal import *
import pysam
from collections import Counter

#######################
def get_consensus_one_allele(data): 
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

######################
numalleles=int(sys.argv[2])
mainfolder_name='/'.join(contig_read_mappings.split('/')[:-2])      #everything but bam file name and sample name
sp_name=contig_read_mappings.split('/')[-2:-1]      #just sample name
final_seqs={}
outfile = open(mainfolder_name+'/loci/'+sp_name+'.fa')

bamfile = pysam.AlignmentFile(sys.argv[1], "rb" )       #open bamfile
for pileupcolumn in samfile.pileup():       #go through each position for each contig
    if pileupcolumn.reference_id not in final_seqs:
        final_seqs[pileupcolumn.reference_id]=[]
    if int(pileupcolumn.n) > 3:     #min 3 reads aligned
        bases_aligned = [pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in pileupcolumn.pileups]     #this is the pileup
        consensus_base = get_consensus(bases_aligned)            
    else:
        consensus_base = 'N'
    final_seqs[pileupcolumn.reference_id].append(consensus_base)

for contig,seq in final_seqs.iteritems():
    outfile.write('>'+contig)
    outfile.write(''.join(seq))

outfile.close()