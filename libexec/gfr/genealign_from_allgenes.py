#! /usr/local/bin/python
import os
import subprocess
import sys
from StringIO import StringIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob

def read_resort(species,seqdict,numalleles):
    print 'Reading data: '+species
    for seq_record in SeqIO.parse(f, "fasta"):
        allele=[]
        gene=list(seq_record.id[:])
        allele.append("".join(species.split('.')[:-1]))
        if numalleles>1:
            allele.append('_')
            allele.append(gene.pop(-2))
            allele.append(gene.pop())
        gene=''.join(gene)
        if gene not in seqdict:
            seqdict[gene]=list()
        seq_record.id=''.join(allele)       #2pspecies
        seqdict[gene].append(seq_record)    
    return seqdict

##########################################################
fafiles = glob.glob("*.fa")
print fafiles
seqdict=dict()
for f in fafiles:
    seqdict = read_resort(f,seqdict,int(sys.argv[1]))        #gene:(SeqRecord(species,sequence))

for gene, records in seqdict.iteritems():
    SeqIO.write(records, gene+'.fa', "fasta")
