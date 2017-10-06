#!/usr/bin/env python2
import Bio
from Bio import SeqIO
import os
import sys
from itertools import tee

inFile  = sys.argv[1]
missing = int(sys.argv[2])
fileName = os.path.basename(inFile)
outName = fileName.replace("_align.fa","_align_RawWithEmptySeqs.fa")
inPath = os.path.dirname(os.path.abspath(inFile))

nono=["N","n","-"]

rawAlign = Bio.SeqIO.parse(inFile, 'fasta')
rawAlign, rawCounter = tee(rawAlign)
rawCount = len(list(rawCounter))

filtered = (rec for rec in rawAlign if any(ch not in nono for ch in rec.seq))
filtered,filterCount = tee(filtered)
newCount = len(list(filterCount))

if rawCount - newCount > missing:
    os.rename(inFile,inPath+"/lostLoci/"+fileName)
else:
    if newCount < rawCount:
        os.rename(inFile,inPath+"/rawAlignmentsWithEmptyTaxa/"+outName)
        Bio.SeqIO.write(filtered, inFile, 'fasta')
