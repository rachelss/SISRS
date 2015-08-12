#!/usr/bin/env python2
import os
import sys
import cPickle
from collections import Counter

#get combined pileup info
def getallbases(path):
    allbases=dict()
    loci={}
    for fi in os.listdir(path):
        if fi.endswith("pileups"):
            filein=open(path+'/'+fi,'r')
            for line in filein:
                splitline=line.split()
                if len(splitline)>4:
                    node,pos,ref,num,bases,qual=line.split()
                    bases=bases.replace('.',ref) #insert ref base
                    bases=bases.replace(',',ref)
                    bases=bases.upper() #everything in uppercase
                    bases=list(bases)
                    
                    loc=node+'/'+pos
                    if loc in allbases:
                        allbases[loc].extend(bases)
                    else:
                        allbases[loc]=bases
                    
                    loci[node]=int(pos)
            filein.close()
    
    return allbases,loci

def remove_extra(base_list):
    bases=['A','C','G','T']
    new_base_list=[]
    ibase_list = iter(base_list)
    for b in ibase_list:
        if b in bases:
            new_base_list.append(b)         #get base
        elif b=='+':                        #skip insertions
            i = int(ibase_list.next())
            while i>0:
                z=ibase_list.next()
                i = i-1
        elif b=='-':                        #get deletion
            i = int(ibase_list.next())
            while i>0:
                z=ibase_list.next()
                i = i-1
            new_base_list.append('-')
        elif b=='^':                        #skip read qual noted at end of read
            z=ibase_list.next()
    
    return new_base_list

###############################################
path=sys.argv[1]

allbases,loci=getallbases(path)      #dictionary of combined pileups - locus/pos:bases(as list)
for pos in allbases:
    bases = remove_extra(allbases[pos])             #remove indel info
    b = Counter(bases).most_common()[0][0]
    if len(b)>0:
        allbases[pos] = b
    else:
        allbases[pos] = 'N'

fasta_dict={l:['N']*n for l,n in loci.iteritems()}
for locus_pos,base in allbases.iteritems():
    locus,pos = locus_pos.split('/')
    fasta_dict[locus][int(pos)-1] = base

output = open(path+'/contigs.fa', 'wb')
for l,seq in fasta_dict.iteritems():
    output.write('>'+l+"\n"+"".join(seq)+"\n")    
output.close()