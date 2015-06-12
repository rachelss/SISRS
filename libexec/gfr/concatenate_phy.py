#!/usr/bin/env python2

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter
import os
import linecache
from collections import defaultdict
import glob

######################

phy_files_file = open(sys.argv[1]+'/phy_files.txt','r')
phy_files = phy_files_file.readlines()
phy_files_file.close()

allspecies = [os.path.basename(sp) for sp in sys.argv[2:]]
cat_data = {s:list() for s in allspecies}

for f in phy_files:
    data = SeqIO.to_dict(SeqIO.parse(f.rstrip(), "phylip-relaxed"))
    for k in allspecies:
        if k in data:
            cat_data[k].extend(list(data[k].seq))
        else:
            cat_data[k].extend(['N'] * len(data[data.keys()[0]].seq))

datalist=[]
for k,v in cat_data.iteritems():
    seq = SeqRecord(Seq(''.join(v)), id=k)
    datalist.append(seq)

SeqIO.write(datalist,sys.argv[1]+'/concat_loci.phylip-relaxed', "phylip-relaxed")


