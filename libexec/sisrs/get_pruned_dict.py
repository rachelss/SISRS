#!/usr/bin/env python2

"""Take multiple mpileup files and determine the base at each position of the reference genome

    arguments:
        path: folder containing mpileup files ending in pileups
        minread: number of reads at a position required to call a base
        thresh: proportion of reads that must be one base for calling to occur

    output:
        path/pruned_dict.pkl : contains pickled dictionary of position:base

"""

import sys
import cPickle
from collections import Counter
import glob
import string
import re

#get combined pileup info
def getallbases(path,minread,thresh):
    assert len(glob.glob1(path,"*.pileups"))==1,'More than one pileup file in'+path
    allbases=dict()
    with open (path+'/'+os.path.basename(path)+'.pileups',"r") as filein:
        for line in filein:
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                cleanBases=getCleanList(ref,bases)  #Get clean bases where * replaced with -
                assert len(cleanBases) == int(num), 'bases are being counted incorrectly: '+ str(bases) + ' should have '+str(num)+' bases, but it is being converted to '+"".join(cleanBases)
                #Extract most common base and its count
                singleBase=(Counter(cleanBases).most_common()[0][0])
                if singleBase == '*':
                    singleBase == '-'
                counts=int((Counter(cleanBases).most_common()[0][1]))

                if counts < minread:
                    finalBase='N'
                else:
                    finalBase=singleBase

                if counts / float(len(bases)) >= thresh:
                    finalBase=singleBase
                else:
                    finalBase='N'
                if finalBase not 'N':
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
            new_base_list.append(b)     #Get base
        elif b in indels:                   #skip indels
            i = int(ibase_list.next())
            j = str(ibase_list.next())
            if str.isdigit(j):
                skip=int(str(i)+j)
                while skip>0:
                    z=ibase_list.next()
                    skip=skip-1
            else:
                while i>1:
                    z=ibase_list.next()
                    i = i-1
        elif b=='^':                        #skip read qual noted at end of read
            z=ibase_list.next()

    return new_base_list

###############################################
if __name__ == "__main__":
    path=sys.argv[1]
    minread=int(sys.argv[2])
    thresh=float(sys.argv[3])

    allbases=getallbases(path,minread,thresh)      #dictionary of combined pileups - locus/pos:bases(as list)
    if len(allbases)==0:
        print 'No data for '+path
        sys.exit(1)      #dictionary of combined pileups - locus/pos:bases(as list)

    output = open(path+'/pruned_dict.pkl', 'wb')
    cPickle.dump(allbases, output, cPickle.HIGHEST_PROTOCOL)
    output.close()
