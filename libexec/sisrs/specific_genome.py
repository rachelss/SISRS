#!/usr/bin/env python2
#This script assumes a single pileup file per taxon
#Add an assert statement to test this
import os
import sys
import cPickle
from collections import Counter
from Bio import SeqIO
import glob

#get combined pileup info
def getallbases(path):
    assert len(glob.glob1(path,"*.pileups"))==1
    allbases=dict()
    with open (path+'/*pileups',"r") as filein:
        for line in filein:
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                bases=getCleanList(ref,bases)
                finalBase=(Counter(bases).most_common()[0][0])
                allbases[loc]=finalBase
    return allbases

def getCleanList(ref,bases):
    bases=bases.replace('.',ref) #insert ref base
    bases=bases.replace(',',ref)
    bases=bases.upper() #everything in uppercase
    bases=list(bases)

    okbases=['A','C','G','T','*']
    indels=['+','-']

    new_base_list=[]
    ibase_list = iter(bases)
    for b in ibase_list:
        if b in okbases:
            if b=='*':
                new_base_list.append('D')   #Replace deletions with Ds as placeholder
            else:
                new_base_list.append(b)     #Get base
        elif b in indels:                   #skip indels
            i = int(ibase_list.next())
            while i>0:
                z=ibase_list.next()
                i = i-1
        elif b=='^':                        #skip read qual noted at end of read
            z=ibase_list.next()

    return new_base_list


###############################################
path=sys.argv[1]
contig_file = sys.argv[2]

allbases=getallbases(path)      #dictionary of combined pileups - locus/pos:bases(as list)
if len(allbases)==0:
    print 'No data for '+path
    sys.exit(1)

# Read contig fasta file into dictionary with sequence ID as the key
contig_handle = open(contig_file, "r")
fasta_seq = SeqIO.parse(contig_handle, 'fasta')
fasta_dict = {read.id:list(str(read.seq)) for read in fasta_seq}
contig_handle.close()

for locus_pos,base in allbases.iteritems():
    locus,pos = locus_pos.split('/')
    fasta_dict[locus][int(pos)-1] = base

output = open(path+'/contigs.fa', 'wb')
for l,seq in fasta_dict.iteritems():
    output.write('>'+str(l)+"\n"+"".join(seq)+"\n")
output.close()
