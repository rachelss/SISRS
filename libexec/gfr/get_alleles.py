#!/usr/bin/env python2
#take fasta file, get pair of alleles
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter
import operator
import math
from decimal import *

#######################        
def factorial_div(numerator, denominator):   #function for when there's a factorial on the top and bottom of a fraction (always assuming the number in the numerator is greater than or equal to the number in the denominator)
    if numerator == denominator:            #shortcut for this case is to cancel out part of the fraction
        return 1                           #e.g. 5! / 3! = 5*4*3*2*1 / 3*2*1 = 5*4
    else:
        return numerator * factorial_div(numerator-1, denominator)

def get_maj_consensus(data):
    pos_bases_counts = Counter(data).most_common()
    c_base = pos_bases_counts[0][0]
    return c_base

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
    
    return c_base, len(data)

def likelihoodtest(listbases,reads1,reads2):  #function to do the likelihood test of Hohenlohe 2010 - takes a sequence in list form
    base_freq=Counter(listbases).most_common(5)
    phasing = 'p'
    if len(reads1) == 0:
        reads1 = [i for i,b in enumerate(list_bases)]
    
    if len(base_freq) == 0:
        return 'n','n',len(listbases),phasing,reads1,reads2
    if len(base_freq) == 1:
        return base_freq[0][0],base_freq[0][0],len(listbases),phasing,reads1,reads2
    else:
        #function from Hohenlohe 2010
        part1a = factorial_div(len(listbases),base_freq[0][1])   #use shortcut method for factorials in fractions
        otherallelecounts = [factorial(c[1]) for i,c in enumerate(base_freq) if i>0]
        part1b = 1/Decimal(reduce(operator.mul, otherallelecounts, 1))
        homopart2=(Decimal(base_freq[0][1])/len(listbases))**(base_freq[0][1])
        homopart3=(Decimal(len(listbases)-base_freq[0][1])/(3*len(listbases)))**(len(listbases)-base_freq[0][1])
        if len(base_freq)>1:
            hetpart2=(Decimal(base_freq[0][1]+base_freq[1][1])/(2*len(listbases)))**(base_freq[0][1]+base_freq[1][1])
        if len(base_freq) == 2:      #anything to the 0 = 1, but the computer needs help with this
            hetpart3=1
        elif len(base_freq) == 3:
            hetpart3=(Decimal(base_freq[2][1])/(2*len(listbases)))**(base_freq[2][1])
        else:
            hetpart3=(Decimal(base_freq[2][1]+base_freq[3][1])/(2*len(listbases)))**(base_freq[2][1]+base_freq[3][1])
        probhomo=part1a*part1b*homopart2*homopart3
        probhet=part1a*part1b*hetpart2*hetpart3
        probhomo+=Decimal(0.000001)      #compensate for limitations of computer
        probhet+=Decimal(0.000001)
        if probhomo=="inf":
            probhomo=Decimal(1.7976931348623157e+308)
        if probhet=="inf":
            probhet=Decimal(1.7976931348623157e+308)
#        print probhet
        ln1=(math.log(probhomo))
        ln2=(math.log(probhet))
        likelihood=Decimal(2*(ln1-ln2))
        if likelihood>=3.84:                 #alpha = 0.05 for test
            return base_freq[0][0], base_freq[0][0],len(listbases),phasing,reads1,reads2             #return correct bases to go into allele
        elif likelihood<=(-3.84):
            newreads1=[i for i,b in enumerate(list_bases) if b is base_freq[0][0]]
            newreads2=[i for i,b in enumerate(list_bases) if b is base_freq[1][0]]
            oneone = set(newreads1).intersection(set(reads1))
            onetwo = set(newreads1).intersection(set(reads2))
            twoone = set(newreads2).intersection(set(reads1))
            twotwo = set(newreads2).intersection(set(reads2))
            if not oneone and not onetwo and not twoone and not twotwo:     #can't link snps
                phasing = 'u'
                return base_freq[0][0], base_freq[1][0],len(listbases),phasing,list(newreads1),list(newreads2)
            elif (len(onetwo)+len(twoone))<(len(oneone)+len(twotwo)):
                phasing = 'p'
                return base_freq[0][0], base_freq[1][0],len(listbases),phasing,list(newreads1),list(newreads2)
            elif (len(onetwo)+len(twoone))>(len(oneone)+len(twotwo)):
                phasing = 'p'
                return base_freq[1][0], base_freq[0][0],len(listbases),phasing,list(newreads2),list(newreads1)
            else:
                phasing = 'u'
                return base_freq[0][0], base_freq[1][0],len(listbases),phasing,list(newreads1),list(newreads2)
        else:
            return 'n','n',len(listbases),phasing,reads1,reads2

######################
bases = ['A','C','G','T','a','c','g','t']
min_bases=1

handle = open(sys.argv[1], "rU")
allseqs = list(SeqIO.parse(handle, "fasta"))
handle.close()

final_seq=[]

for i in range(len(allseqs[0])):
    data = [seq[i] for seq in allseqs if seq[i] in bases]
    if len(data)>=min_bases:
        b = get_consensus(data)
        final_seq.append(b)
    elif 0<len(data)<min_bases:
        final_seq.append('N')
        


outfile = open(sys.argv[1].replace('_align.fa','_alleles.fa'),'w')
outfile.write('>'+sys.argv[1].replace('_align.fa',''))
outfile.write("".join(final_seq))
outfile.close()