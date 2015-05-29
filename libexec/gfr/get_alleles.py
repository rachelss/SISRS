#!/usr/bin/env python2
#take fasta file, get pair of alleles
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter

#######################
def get_consensus(data):
    pos_bases_counts = Counter(data).most_common()
    c_base = pos_bases_counts[0][0]
    return c_base

######################
bases = ['A','C','G','T','a','c','g','t']

handle = open(sys.argv[1], "rU")
allseqs = list(SeqIO.parse(handle, "fasta"))
handle.close()

final_seq=[]

for i in range(len(allseqs[0])):
    data = [seq[i] for seq in allseqs if seq[i] in bases]
    if len(data)>0:
        b = get_consensus(data)
        final_seq.append(b)

outfile = open(sys.argv[1].replace('_align.fa','_alleles.fa'),'w')
outfile.write('>'+sys.argv[1].replace('_align.fa',''))
outfile.write("".join(final_seq))
outfile.close()