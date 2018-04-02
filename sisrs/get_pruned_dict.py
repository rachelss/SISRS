#!/usr/bin/env python3

"""Take multiple mpileup files and determine the base at each position of the reference genome

    arguments:
        path: folder containing mpileup files ending in pileups
        contig_dir: folder containing assembled composite genome
        minread: number of reads at a position required to call a base
        thresh: proportion of reads that must be one base for calling to occur

    output:
        path/<SP>_LocList: list of bases for each position
"""

import sys
from collections import Counter
import glob
import string
import re
import os
from .specific_genome import getCleanList
from collections import defaultdict

#get combined pileup info
def getallbases(path,contig_dir,minread,thresh):
    basePath=os.path.dirname(path)
    assert len(glob.glob1(path,"*.pileups"))==1,'More than one pileup file in'+path
    speciesDict=defaultdict(lambda: 'N')
    minPenalty=0
    threshPenalty=0
    bothPenalty=0
    with open (path+'/'+os.path.basename(path)+'.pileups',"r") as filein:
        for line in iter(filein):
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                cleanBases=getCleanList(ref,bases)  #Get clean bases where * replaced with -
                finalBase,minPenalty,threshPenalty,bothPenalty=getFinalBase_Pruned(cleanBases,minread,thresh,minPenalty,threshPenalty,bothPenalty)
                speciesDict[loc] = finalBase

    printSpecies = open(path+"/"+os.path.basename(path)+'_LocList', 'w')
    with open(contig_dir+"/contigs_LocList") as f:
        for line in f:
            printSpecies.write(speciesDict[line.strip()]+'\n')
    f.close()
    printSpecies.close()

    c = Counter(speciesDict.values())
    nCount = c['N']
    siteCount = len(speciesDict) - nCount
    sitePercent = format((float(siteCount)/len(speciesDict))*100,'.2f')
    nPercent = format((float(nCount)/len(speciesDict))*100,'.2f')
    print("Of "+ str(len(speciesDict)) + " positions, " + os.path.basename(path) + " has good calls for " + str(siteCount) + " sites (" + sitePercent +"%). There were " + str(nCount) + " N calls ("+ nPercent + "%).",flush=True)
    print("Of " + str(nCount) + " Ns, " + os.path.basename(path) + " lost " + str(threshPenalty) + " via homozygosity threshold, " + str(minPenalty)  +" from low coverage, and " + str(bothPenalty) + " from both. "+ str(nCount - threshPenalty - minPenalty - bothPenalty) + " sites had no pileup data.\n",flush=True)

    return siteCount

def getFinalBase_Pruned(cleanBases,minread,thresh,minPenalty,threshPenalty,bothPenalty):
    singleBase=(Counter(cleanBases).most_common()[0][0])
    if singleBase == '*':
        singleBase = '-'
    counts=int((Counter(cleanBases).most_common()[0][1]))

    if counts >= minread and counts/float(len(cleanBases)) >= thresh:
        finalBase=singleBase
    else:
        finalBase='N'
        if counts < minread and counts/float(len(cleanBases)) < thresh:
            bothPenalty+=1
        elif counts < minread:
                minPenalty+=1
        elif counts/float(len(cleanBases)) < thresh:
                threshPenalty+=1

    return finalBase,minPenalty,threshPenalty,bothPenalty
###############################################
def main(path, contig_dir, minread, thresh):

    allbases=getallbases(path,contig_dir,minread,thresh)      #dictionary of combined pileups - locus/pos:bases(as list)
    if allbases==0:
        print('No data for '+path,flush=True)
        sys.exit(1)

if __name__ == '__main__':
    path = sys.argv[1]
    contig_dir = sys.argv[2]
    minread = int(sys.argv[3])
    thresh = float(sys.argv[4])
    main(path, contig_dir, minread, thresh)
