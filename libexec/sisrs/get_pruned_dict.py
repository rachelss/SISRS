#!/usr/bin/env python2
import sys
import cPickle
from collections import Counter
import glob
import string
import re

#get combined pileup info
def getallbases(path):
    allbases=dict()
    for fi in glob.glob(path+'/*pileups'):
        filein=open(fi,'r')
        for line in filein:
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                bases=bases.replace('.',ref) #insert ref base
                bases=bases.replace(',',ref)
                bases=bases.replace('*','D')    #using D is easier than *
                bases=bases.upper() #everything in uppercase
                bases = remove_extra(bases)

                assert len(bases) == num, 'bases are being counted incorrectly '+\
                                          "".join(bases)+' should only be '+str(num)+' bases'

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
        counts = Counter(bases).most_common(1)
        if counts[0][1] / float(len(bases)) >= thresh:
            base=counts[0][0]
        else:
            base='N'
                
    return base

def remove_extra(base_str):
    bases=['A','C','G','T','D']     #D indicates del
    new_base_list=[]
    ibase_list = iter(re.findall('[A-Z]+|\+[0-9]+|\-[0-9]+|\^.',base_str))  #group by bases, indels, read start (need this to get rid of following ascii code), don't worry about $
    for b in ibase_list:
        if b[0] in string.ascii_uppercase:     #allows for any other characters - filter out later (in case of ZAAA)
            new_base_list.extend(b)         #get bases
        elif b[0]=='+' or b[0]=='-':        #skip indels
            i = int(b[1:])                  #account for ints >=10
            b2=ibase_list.next()
            if len(b2) > i:
                new_base_list.extend(b2[i:])    #can have real bases after indel bases
        elif b[0]=='^':                        #skip read qual noted at end of read
            continue
    new_base_list = [i for i in new_base_list if i in bases]    #filter unknown bases
    
    return new_base_list        #list of individual bases

def test_remove_extra(teststring,outstring):
    o = remove_extra(teststring)
    assert o == list(outstring), 'Test failure: result is '+"".join(o)+' not '+outstring
    
###############################################
path=sys.argv[1]
minread=int(sys.argv[2])
thresh=float(sys.argv[3])

test_remove_extra('A+1CAAA-10*******AAA^--1*AA-2CAA-2**A-2***','AAAAAAAAD')

allbases=getallbases(path)      #dictionary of combined pileups - locus/pos:bases(as list)
for pos in allbases:
    base = determine_base(bases,minread,thresh)     #determine if sufficient data and threshold met for calling allele
    if base is 'D':
        base='-'
    allbases[pos] = base

output = open(path+'/pruned_dict.pkl', 'wb')
cPickle.dump(allbases, output, cPickle.HIGHEST_PROTOCOL)
output.close()