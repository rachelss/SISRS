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
import itertools
from itertools import izip_longest

#get combined pileup info
def getallbases(posList,minread,thresh):
    assert len(glob.glob1(path,"*.pileups"))==1,'More than one pileup file in'+path

    keys = posList
    values = ['N'] * len(posList)
    speciesDict = dict(zip(keys, values))

    with open (path+'/'+os.path.basename(path)+'.pileups',"r") as filein:
        for line in iter(filein):
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                cleanBases=getCleanList(ref,bases)  #Get clean bases where * replaced with -
                finalBase=getFinalBase_Pruned(cleanBases,minread,thresh)
                speciesDict[loc] = finalBase
    printSpecies = open(path+"/"+os.path.basename(path)+'_LocList', 'w')
    for item in posList:
        print>>printSpecies, speciesDict[item]
    printSpecies.close()

    valueList=[]
    for key,value in speciesDict.iteritems():
        valueList.append(value)
    nCount = valueList.count("N")
    siteCount = len(speciesDict) - nCount
    sitePercent = format((float(siteCount)/len(speciesDict))*100,'.2f')
    nPercent = format((float(nCount)/len(speciesDict))*100,'.2f')
    sys.stdout.write("Of "+ str(len(speciesDict)) + " positions, " + os.path.basename(path) + " has good calls for " + str(siteCount) + " sites (" + sitePercent +"%). There were " + str(nCount) + " N calls ("+ nPercent + "%)\n")
    return siteCount

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
    posList=[]
    with open("/data3/schwartzlab/bob/Mammal_SISRS/SISRS_Runs/A/premadeoutput/contigs_LocList") as f:
        for line in f:
            posList.append(line.strip())
    sys.stdout.write("List read. Size: "+str(sys.getsizeof(posList)))

    #Generate species-specific posList
    siteCount=getallbases(posList,minread,thresh)      #dictionary of combined pileups - locus/pos:bases(as list)
    if siteCount == 0:
        print 'No data for '+path
        sys.exit(1)
