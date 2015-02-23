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
                print bases, len(bases),str(Counter(bases).most_common(1)[0][1] / float(len(bases)))
                if(len(bases)< minread):                                                #not enough info
                    allbases[str(bamfile.getrname(pileupcolumn.tid))+str(pileupcolumn.pos)]='N'      #node_pos:'N'
                    print '1: N'
                elif Counter(bases).most_common(1)[0][1] / float(len(bases)) >= thresh: #enough of one
                    allbases[str(bamfile.getrname(pileupcolumn.tid))+'_'+str(pileupcolumn.pos)]=Counter(bases).most_common(1)[0][0] #node_pos:base
                    print '2: ',Counter(bases).most_common(1)[0][0]
                else:       #het or lots of error
                    allbases[str(bamfile.getrname(pileupcolumn.tid))+str(pileupcolumn.pos)]='N'      #node_pos:'N'
                    print '1: N'
        
            bamfile.close()
    
    return allbases

###############################################
allbases=getallbases(sys.argv[1],sys.argv[2],sys.argv[3])      #dictionary of combined pileups - locus/pos:bases(as list)
output = open(sys.argv[1]+'/pruned_dict.pkl', 'wb')
cPickle.dump(allbases, output, cPickle.HIGHEST_PROTOCOL)
output.close()