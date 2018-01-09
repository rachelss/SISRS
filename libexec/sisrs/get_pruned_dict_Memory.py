#!/usr/bin/env python2

"""Take multiple mpileup files and determine the base at each position of the reference genome

    arguments:
        path: folder containing mpileup files ending in pileups
        assembler: Assembler used for SISRS or 'premade'
        minread: number of reads at a position required to call a base
        thresh: proportion of reads that must be one base for calling to occur
"""

import sys
from collections import Counter
import glob
import string
import re
import os
from specific_genome import getCleanList

#get combined pileup info
def getallbases(posList,minread,thresh):
    assert len(glob.glob1(path,"*.pileups"))==1,'More than one pileup file in'+path
    speciesList = ['N'] * len(posList)
    with open (path+'/'+os.path.basename(path)+'.pileups',"r") as filein:
        for line in filein:
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                pos=siteList.index(loc)
                cleanBases=getCleanList(ref,bases)  #Get clean bases where * replaced with -
                finalBase=getFinalBase_Pruned(cleanBases,minread,thresh)
                speciesList.insert(pos,finalBase)
    printSpecies = open(path+"/"+os.basename(path)+'_LocList', 'w')
    for item in speciesList:
        print>>printSpecies, item
    printSpecies.close()
    nCount = len(posList) - speciesList.count("N")
    siteCount = len(posList) - nCount
    print "Of "+ len(posList) + "positions, " + os.basename(path) + "has non-N sites for " + str(siteCount) + "sites."
    return nCount

def getFinalBase_Pruned(cleanBases,minread,thresh):
    singleBase=(Counter(cleanBases).most_common()[0][0])
    if singleBase == '*':
        singleBase = '-'
    counts=int((Counter(cleanBases).most_common()[0][1]))

    if counts >= minread and counts/float(len(cleanBases)) >= thresh:
        finalBase=singleBase
    else:
        finalBase='N'

    return finalBase
###############################################
if __name__ == "__main__":

    #Read in arguments
    path=sys.argv[1]
    basePath=os.path.dirname(path)
    assembler=sys.argv[2]
    minread=int(sys.argv[3])
    thresh=float(sys.argv[4])

    #Read in PosList
    with open(basePath+"/"+assembler+"output/contigs_PosList") as f:
        posList = f.read().splitlines()

    #Generate species-specific posList
    allbases=getallbases(posList,minread,thresh)      #dictionary of combined pileups - locus/pos:bases(as list)
    if allbases == 0:
        print 'No data for '+path
        sys.exit(1)
