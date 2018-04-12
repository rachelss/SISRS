#!/usr/bin/env python2

""" Subsample fastq data using reservoir sampling
    assumes reads are 100 bp
    samples read pairs or unpaired
    
    arguments:
        number of pairs of reads (or half the number of total reads)
        folder in which to find fastq
        label for forward reads
        label for reverse reads
    
    output:
        file with paired samples
        file with unpaired samples
            in subsamples folder

    run: parallel --jobs "${PROCESSORS}" "python ${DIR}/libexec/sisrs/sub_sample_for_velvet_unshuff.py ${LEFTREADS} {}" ::: "${FOLDERLISTA[@]}"
"""

import sys
import random
import glob
import os

def main(taxonfolder, left_reads, fileid, fileid2):

    N = int(left_reads)
    #N = int(sys.argv[1])     #number of left reads required per sample for 10x total coverage; (10*genome size) / (read size*2*num_samples)
    samplep,sampleu = [],[]
    i=0

    #taxonfolder=sys.argv[2]
    mainfolder=os.path.dirname(os.path.abspath(taxonfolder))
    taxon_files = glob.glob(taxonfolder+'/*fastq')
    taxon=os.path.basename(taxonfolder)

    #fileid = sys.argv[3]
    #fileid2 = sys.argv[4]

    paired1 = [f for f in taxon_files if fileid in f]
    print('paired files: '+" ".join(paired1))
    unpaired = [f for f in taxon_files if fileid not in f and fileid2 not in f]
    print('unpaired files: '+" ".join(unpaired))

    for fq in paired1:
        f1=open(fq, 'r')
        f2=open(fq.replace(fileid,fileid2), 'r')
        
        for line in f1:
            if not line: break      #have less than desired coverage
            read_info=[line,next(f1),next(f1),next(f1),next(f2),next(f2),next(f2),next(f2)]
            if i < N:        
                samplep.append(read_info)
            elif random.random() < N/float(i+1):
                replace = random.randint(0,N-1)
                samplep[replace] = read_info
            i+=1
        f1.close()
        f2.close()

    for fq in unpaired:
        f1=open(fq, 'r')
        for line in f1:
            if not line: break      #have less than desired coverage
            read_info=[line,next(f1),next(f1),next(f1)]
            try:
                read_info.extend([next(f1),next(f1),next(f1),next(f1)])
            except:
                break
            if i < N:        
                sampleu.append(read_info)
            elif random.random() < N/float(i+1):   
                replace = random.randint(0,N-1)
                if replace<len(samplep):
                    del samplep[replace]
                else:
                    del sampleu[replace-len(samplep)]
                sampleu.append(read_info)
            i+=1
        f1.close()

    if len(samplep)>0:            
        f = open(mainfolder+'/subsamples/'+taxon+'subsampledp.fastq','w')
        for read_info in samplep:
            f.write("".join(read_info))
        f.close()
        print(str(len(samplep))+' paired reads sampled')

    if len(sampleu)>0:
        f = open(mainfolder+'/subsamples/'+taxon+'subsampledu.fastq','w')
        for read_info in sampleu:
            f.write("".join(read_info))
        f.close()
        print(str(len(sampleu)*2)+' unpaired reads sampled')
