#!/usr/bin/env python2
#take fasta file, get pair of alleles
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter,defaultdict
import operator
import math
from decimal import *
import os.path

#######################        
def factorial_div(numerator, denominator):   #function for when there's a factorial on the top and bottom of a fraction (always assuming the number in the numerator is greater than or equal to the number in the denominator)
    if numerator == denominator:            #shortcut for this case is to cancel out part of the fraction
        return 1                           #e.g. 5! / 3! = 5*4*3*2*1 / 3*2*1 = 5*4
    else:
        return numerator * factorial_div(numerator-1, denominator)

def trimends(allele):
    if allele[0]=='N':
        while allele[0]=='N':
            allele.pop(0)
    if allele[-1]=='N':
        while allele[-1]=='N':
            allele.pop()
    return allele

def get_maj_consensus(data):
    pos_bases_counts = Counter(data).most_common()
    c_base = pos_bases_counts[0][0]
    return c_base

def get_80_20_consensus(data):      #takes list of bases, returns consensus based on 80/20 rule
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

def get_one_allele(basereads,min_bases):      #takes list of info for each position when info is dict(base:[reads], base:[reads]), returns an allele
    final_allele=[]
    for pos in basereads:
        baselist = [[b]*len(ids) for b,ids in pos.iteritems()]
        baselist = [item for inner_list in baselist for item in inner_list]
        base = get_80_20_consensus(baselist)
        if len(baselist)>=min_bases:
            final_allele.append(base)
        else:
            final_allele.append('N')
            
    final_allele=trimends(final_allele)
            
    return("".join(final_allele))
    
def likelihoodtest(bpreads):        #returns one base (one or two copies) or two bases
    bpfreq = {b: float(len(r)) for b,r in bpreads.iteritems()}
    bpfreqlist = sorted(bpfreq.items(), key=operator.itemgetter(1),reverse=True)      #sort bases by frequency and put in list
    n=[i[1] for i in bpfreqlist]
    
    A=n[0]/sum(n)
    A=math.log(A)
    A=n[0]*A
    Bn=sum(n[1:])
    Bd=3*sum(n)
    B=Bn*math.log(Bn/Bd)
    Cn=sum(n[0:2])
    Cd=2*sum(n)
    C=Cn*math.log(Cn/Cd)
    if len(n)>2:
        Dn=sum(n[2:])
        Dd=2*sum(n)
        D=Dn*math.log(Dn/Dd)
    else:
        D=0
    LR = 2*(A+B-C-D)
    
    if LR>=3.84:                 #alpha = 0.05 for test
        return bpfreqlist[0][0], bpfreqlist[0][0]             #return correct bases to go into allele
    elif LR<=(-3.84):
        return bpfreqlist[0][0], bpfreqlist[1][0]
    else:
        return bpfreqlist[0][0],'n'
        
def phase(bpreads,pastreads1,pastreads2):   #returns two bases in correct order given current and previous snp reads
    bpfreq = {b: float(len(r)) for b,r in bpreads.iteritems()}
    bpfreqlist = sorted(bpfreq.items(), key=operator.itemgetter(1),reverse=True)
    bases = [bpfreqlist[0][0],bpfreqlist[1][0]]
    reads1=bpreads[bases[0]]
    reads2=bpreads[bases[1]]

    d = set(reads1).intersection(set(pastreads1))
    e = set(reads2).intersection(set(pastreads2))
    f = set(reads1).intersection(set(pastreads2))
    g = set(reads2).intersection(set(pastreads1))

    if (len(d)+len(e))<(len(f)+len(g)):          #if there are more shared if the alleles were swapped, swap
        return bases[1],bases[0],reads2,reads1
    
    else:     #if no linkage, uncertain linkage, or correct linkage leave as is
        return bases[0],bases[1],reads1,reads2

def get_two_alleles(basesreads,min_bases):    #returns two strings = alleles
    allele1,allele2=[],[]
    reads1,reads2=[],[]
    for pos in basereads:
        baselist = [[b]*len(ids) for b,ids in pos.iteritems()]
        baselist = [item for inner_list in baselist for item in inner_list]
        
        if len(baselist)>=min_bases:        
            if len(pos)==1:
                allele1.append(pos.keys()[0])
                allele2.append(pos.keys()[0])
            else:
                pos_allele1,pos_allele2 = likelihoodtest(pos)
                if (pos_allele2 != pos_allele1) and (pos_allele2 != 'n'):
                    pos_allele1,pos_allele2, reads1, reads2 = phase(pos,reads1,reads2)
                allele1.append(pos_allele1)
                allele2.append(pos_allele2)
        else:
            allele1.append('N')
            allele2.append('N')
            
    allele1 = trimends(allele1)
    allele2 = trimends(allele2)
    
    return "".join(allele1),"".join(allele2)

def base_reads(seqlist):
    allinfo=[]
    bases = ['A','C','G','T','a','c','g','t']
    for i in range(len(seqlist[0])):
        reads = defaultdict(list)
        for seq in seqlist:
            if seq[i] in bases:
                reads[seq[i]].append(seq.id)
        allinfo.append(reads)
    return allinfo         #[{base:[id,id...],base:[id,id...]}, ]   


######################
if __name__ == '__main__':

    min_bases=int(sys.argv[3])
    sp=sys.argv[1].split('/')[-2].replace('_loci','')        #species = folder name
    num_alleles=int(sys.argv[2])
    
    if os.path.exists(sys.argv[1].replace('.fa','_align.fasta')):
        handle = open(sys.argv[1].replace('.fa','_align.fasta'), "rU")
        allseqs = list(SeqIO.parse(handle, "fasta"))
        handle.close()
    
        basereads = base_reads(allseqs)      #[{base:[id,id...],base:[id,id...]}, ]
        
        if num_alleles==1:
            allele1 = get_one_allele(basereads,min_bases)
        elif num_alleles==2:
            allele1, allele2 = get_two_alleles(basereads,min_bases)
        else:
            print 'You must have either one or two alleles'
            sys.exit(1)
    
    else:     # didn't align b/c only one sequence
        if min_bases==1:
            handle = open(sys.argv[1], "rU")     #read that one seq
            seq=list(SeqIO.parse(handle, "fasta"))
            handle.close()
            
            allele1=str(seq[0].seq)
            num_alleles=1
        else:
            'There was only one read for '+sys.argv[1]
            sys.exit(1)
    
    outfile = open(sys.argv[1].replace('.fa','_alleles.fa'),'w')
    outfile.write('>'+sp+"\n")
    outfile.write(allele1+"\n")
    if num_alleles==2:
        outfile.write('>'+sp+'_2'+"\n")
        outfile.write(allele2+"\n")        
    outfile.close()