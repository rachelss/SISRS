#! /usr/local/bin/python
#take sam file aligned to single contig and produce whole alignment
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
def factorial(n):        #function to do factorials
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)
         
def factorial_div(numerator, denominator):   #function for when there's a factorial on the top and bottom of a fraction (always assuming the number in the numerator is greater than or equal to the number in the denominator)
    if numerator == denominator:            #shortcut for this case is to cancel out part of the fraction
        return 1                           #e.g. 5! / 3! = 5*4*3*2*1 / 3*2*1 = 5*4
    else:
        return numerator * factorial_div(numerator-1, denominator)


def sep_cigar(cigar):
    d=list()
    code=list()
    nums=list()
    for i,c in enumerate(cigar):
        if c.isdigit():
            d.append(c)
        else:
            code.append(c)
            d="".join(d)
            nums.append(int(d))
            d=list()

    return nums, code


def adjust_seq(seq,qual,cigar):
    seq = list(seq)
    qual = list(qual)
    newseq,newqual = [],[]
    nums,code = sep_cigar(cigar)
    for j,num in enumerate(nums):
        if code[j] == 'S':
            if j == 0:
                for z in range(num):
                    b = seq.pop(0)
                    c = qual.pop(0)
            else:
                for z in range(num):
                    b = seq.pop(0)
                    c = qual.pop(0)
                    newseq.append(b)
                    newqual.append(c)
        elif code[j] == 'I':
            for z in range(num):
                b = seq.pop(0)
                c = qual.pop(0)
        elif code[j] == 'D':
            for z in range(num):
                newseq.append('-')
                newqual.append('-')
        elif code[j] == 'M':
            for z in range(num):
                b = seq.pop(0)
                c = qual.pop(0)
                newseq.append(b)
                newqual.append(c)
    
    return newseq,newqual

def likelihoodtest(list_bases,reads1,reads2):  #function to do the likelihood test of Hohenlohe 2010 - takes a sequence in list form
    listbases = [b for b in list_bases if b in ['A','C','G','T','a','g','c','t']]
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
######################
contig_read_mappings=sys.argv[1]
numalleles=int(sys.argv[2])
print contig_read_mappings,
folder_name=contig_read_mappings.split('/')[0]
file_name=contig_read_mappings.split('/')[1]
node_name=file_name.split('.')[0]

totalseq,totalqual=[],[]
#get alignment
samfile=open(contig_read_mappings,'r')    #open file
for line in samfile:      #go through file
    if line.startswith('@'):
        continue
    elif len(line.strip())>0:
        splitline=line.split()
        flag = bin(int(splitline[1]))
        flagl = flag.split('b')
        if len(flagl[1])>=3:
            flag = flagl[1][-3]     #check for mapping
            if flag == '1':
                continue
        if len(flagl[1])>=5:
            flag = flagl[1][-5]     #check for reverse mapping
        else:
            flag = '0'
            
        cigar = splitline[5]
        pos = int(splitline[3])
        seq=splitline[9]
        qual = splitline[10]
            
        seq,qual = adjust_seq(seq,qual,cigar)
        
        while pos>1:
            seq.insert(0, '-')
            qual.insert(0, '-')
            pos=pos-1
        totalseq.append("".join(seq))
        totalqual.append("".join(qual))

numreads = len(totalseq)
if numreads == 0:
    print sys.argv[1]
maxlen=0
for i,j in enumerate(totalseq):
    if len(j)>maxlen:
        maxlen = len(j)

for i,j in enumerate(totalseq):
    plus = maxlen - len(j)
    totalseq[i] = j + ('-' * plus)          #reads have - added to make them all the same length
    totalqual[i] = j + ('-' * plus)

if numalleles == 1:
    print 'one allele'
    allele1,site_counts = [],[]
    for i in range(len(totalseq[0])):
        base1,num = get_consensus([seq[i] for s,seq in enumerate(totalseq) if seq[i] is not '-' and ord(totalqual[s][i]) >39])  #7+32
        allele1.append(base1)
    allele1_seq_r = SeqRecord(Seq(''.join(allele1)), id=node_name)
    alleles = []
    alleles.append(allele1_seq_r)
    SeqIO.write(alleles,folder_name+'/'+node_name+'.fa', "fasta")
else:
    allele1,allele2,site_counts,reads1,reads2 = [],[],[],[],[]
    phasingp='p'
    for i in range(len(totalseq[0])):
        base1, base2,num,phasing,reads1,reads2 = likelihoodtest([seq[i] for s,seq in enumerate(totalseq) if ord(totalqual[s][i]) >39],reads1,reads2)
        allele1.append(base1)
        allele2.append(base2)
        site_counts.append(str(num))
        if phasing == 'u':
            phasingp = 'u'
    allele1_seq_r = SeqRecord(Seq(''.join(allele1)), id=node_name+phasingp+'1')
    allele2_seq_r = SeqRecord(Seq(''.join(allele2)), id=node_name+phasingp+'2')
    alleles = []
    alleles.append(allele1_seq_r)
    alleles.append(allele2_seq_r)
    SeqIO.write(alleles,folder_name+'/'+node_name+'.fa', "fasta")

outfile=open(folder_name+'/'+node_name+'.counts','w')
outfile.write(' '.join(site_counts))
outfile.close()             

#print numreads