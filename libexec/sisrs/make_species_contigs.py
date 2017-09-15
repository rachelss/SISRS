#! /usr/local/bin/python
#take sam file aligned to single contig and produce whole alignment
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter
from os import path

######################
bases = ['A','C','G','T','a','c','g','t']

data = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
contig_names = data.keys()

for k in data:
    data[k]=list(data[k].seq)

pileup = open(sys.argv[2],'r')
for line in pileup:
    splitline = line.split()
    node = splitline[0]
    pos = int(splitline[1])-1
    if len(splitline)>4:
        b=splitline[4].replace('.',splitline[2])
        site = [i for i in b if i in bases]
        if len(site)>0:
            sitecounts = Counter(site).most_common(1)
            data[node][pos] = sitecounts[0][0]
    
datalist=[]
for k,v in data.iteritems():
    seq = SeqRecord(Seq(''.join(v)), id=k)
    datalist.append(seq)
    
SeqIO.write(datalist,path.dirname(sys.argv[2])+'/contigs.fa', "fasta")