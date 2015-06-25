#!/usr/bin/env python2
#takes bam file aligned to ref_genes for one species
#outputs .fa file for each locus for one species

import sys
from decimal import *
from collections import defaultdict
import string
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re

class splitseq:
    def __init__(self,qname,h,a,t):
        self.qname=qname
        self.h=h
        self.a=a
        self.t=t

#################################
infile = sys.argv[1]
mainfolder_name='/'.join(infile.split('/')[:-2])
folder_name='/'.join(infile.split('/')[:-1])      #everything but sam file name

#get contig lengths
contigs = SeqIO.to_dict(SeqIO.parse(mainfolder_name+'/ref_genes.fa', "fasta"))
for k,v in contigs.iteritems():
    contigs[k]=len(v)

reads=defaultdict(list)
samfile = open(infile,'r')
qnames=[]
maxh=0
maxt=0
for line in samfile:
    splitline = line.split()
    qname=splitline[0]
    contig=splitline[2]
    pos=int(splitline[3])
    mapq=int(splitline[4])
    cigar=splitline[5]
    seq=list(splitline[9])
    if mapq>10:
        nums=re.findall('\d+', cigar)
        w=re.findall('\D+', cigar)
        if 'S' == w[0]:
            leftclip=int(nums[0])
        else:
            leftclip=0
        if 'S' == w[-1]:
            rightclip=int(nums[-1])
        else:
            rightclip=0
        
        #get aligned region
        if rightclip==0:
            a=seq[leftclip:]
        else:
            a=seq[leftclip:-rightclip]
        a=['N']*(pos-1)+a
        if pos-1+len(seq) < contigs[contig]:
            a=a+['N']*(contigs[contig]-len(seq)-pos+1)
        
        #get left (head)
        if pos==1:
            h=seq[:leftclip]
            if len(h)>maxh:
                maxh=len(h)
        else:
            h=[]
        
        #get right (tail)
        if pos-1+len(seq)>contigs[contig]:
            t=seq[-rightclip:]
            if len(t)>maxt:
                maxt=len(t)
        else:
            t=[]
        
        #put in dict
        if qname in qnames:
            qname=qname+'_2'
        reads[contig].append(splitseq(qname,h,a,t))    #note: read is already rev-comp if necessary
        qnames.append(qname)

samfile.close()

#add N's to each part of alignment to even things out
for contig,readlist in reads.iteritems():
    outfile = open(folder_name+'/'+contig+'.fa','w')
    for i,r in enumerate(readlist):
        if len(r.h) < maxh:
            reads[contig][i].h = ['N']*(maxh-len(r.h))+r.h
        if len(r.t) < maxt:
            reads[contig][i].t = r.t+['N']*(maxt-len(r.t))   
        outfile.write('>'+r.qname+"\n"+"".join(r.h+r.a+r.t)+"\n")
    outfile.close()