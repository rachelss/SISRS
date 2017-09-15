#!/usr/bin/env python2
from Bio import SeqIO
import os
import sys

path=sys.argv[1]
contigFile=(path+'/contigs.fa')

file = open(path+'/contigs_SeqLength.tsv', "w")

for seq_record in SeqIO.parse(contigFile,"fasta"):
	file.write(str(seq_record.id)+"\t"+str(len(seq_record))+"\n")
file.close()
