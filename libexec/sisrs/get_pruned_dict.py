#!/usr/bin/env python2
import os
import sys
import cPickle
from collections import Counter

#get combined pileup info
def getallbases(path):
    allbases=dict()
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
            filein.close()
    
    return allbases

#prune by whether there's enough info and species are fixed
def determine_base(bases,minread,thresh):
    if(len(bases)< minread): #not enough info
        base='N'
    else:
        counts = Counter(bases)
        if len(counts)==1:
            base=bases[0]
        else:
            if counts.most_common(1)[0][1] / float(len(bases)) >= thresh:
                base=counts.most_common(1)[0][0]
            else:
                base='N'
                
    return base

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
minread=int(sys.argv[2])
thresh=float(sys.argv[3])

allbases=getallbases(path)      #dictionary of combined pileups - locus/pos:bases(as list)
for pos in allbases:
    bases = remove_extra(allbases[pos])             #remove indel info
    base = determine_base(bases,minread,thresh)     #determine if sufficient data and threshold met for calling allele
    allbases[pos] = base

output = open(path+'/pruned_dict.pkl', 'wb')
cPickle.dump(allbases, output, cPickle.HIGHEST_PROTOCOL)
output.close()