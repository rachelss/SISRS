#!/usr/bin/env python2
#takes nexus file, outputs phylip file
#only X species allowed to be missing per site

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter
from os import path
import linecache

######################
bases = ['A','C','G','T','a','c','g','t']
formats = {'nex':'nexus', 'phy':'phylip-relaxed','fa':'fasta'}
fformat = formats[sys.argv[1].split('.')[-1]]
missing=int(sys.argv[2])
data = SeqIO.to_dict(SeqIO.parse(sys.argv[1], fformat))
locline = linecache.getline(sys.argv[1], 8)
locs = locline.split()
species = data.keys()
minsp=len(species)-missing
newlocs=list()

for k in data:
    data[k]=list(data[k].seq)

newdata = {sp: list() for sp in species}
for i in range(len(data[species[0]])):
    site = [data[sp][i] for sp in species if data[sp][i] in bases]
    if len(set(site))>1 and len(site)>=minsp:
        for sp in species:
            newdata[sp].append(data[sp][i])
        newlocs.append(locs[i+1])

   
datalist=[]
for k,v in newdata.iteritems():
    seq = SeqRecord(Seq(''.join(v)), id=k)
    datalist.append(seq)
    
SeqIO.write(datalist,"".join(sys.argv[1].split('.')[:-1])+'_m'+sys.argv[2]+'.phylip-relaxed', "phylip-relaxed")
locfile = open(path.basename(sys.argv[1])+'/locs_m'+sys.argv[2]+'.txt','w')
locfile.write(" ".join(newlocs))
locfile.close()