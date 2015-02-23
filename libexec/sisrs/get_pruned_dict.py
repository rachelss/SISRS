#!/usr/bin/env python2
import os
import sys
import cPickle
from collections import Counter
import pysam

#get combined pileup info
def getallbases(path,minread,thresh):
    allbases=dict()
    for fi in os.listdir(path):
        if fi.endswith("bam"):
            bamfile = pysam.AlignmentFile(path+'/'+fi, "rb" )
            for pileupcolumn in bamfile.pileup():           #only doing one at a time
                basesall=[pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in pileupcolumn.pileups]       #get bases per site
                bases=[b for b in basesall if b in ['a','c','g','t','A','C','G','T']]
                nodepos=str(bamfile.getrname(pileupcolumn.tid))+'_'+str(pileupcolumn.pos)
                if nodepos in allbases:
                    allbases[nodepos].extend(bases)
                else:
                    allbases[nodepos]=bases
            bamfile.close()
            
    for k,v in allbases.iteritems():
        if(len(v)< minread):                                                #not enough info
            allbases[k]='N'      #node_pos:'N'
        elif Counter(v).most_common(1)[0][1] / float(len(v)) >= thresh: #enough of one
            allbases[k]=Counter(v).most_common(1)[0][0] #node_pos:base
        else:       #het or lots of error
            allbases[k]='N'      #node_pos:'N'        
            
    return allbases

###############################################
allbases=getallbases(sys.argv[1],int(sys.argv[2]),float(sys.argv[3]))      #dictionary of combined pileups - locus/pos:bases(as list)
output = open(sys.argv[1]+'/pruned_dict.pkl', 'wb')
cPickle.dump(allbases, output, cPickle.HIGHEST_PROTOCOL)
output.close()