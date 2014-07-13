#!/usr/bin/python
import sys
import random
import glob

N = int(sys.argv[1])     #number of left reads required per sample for 10x total coverage; (10*genome size) / (read size*2*num_samples)
sample = [];
i=0

for fq in glob.glob(sys.argv[2]+"/*shuffled*fastq"):
    print 'Subsampling '+fq
    f = open(fq,'r')
    for line in f:
        if not line: break      #have less than desired coverage
        read_info=[line,f.next(),f.next(),f.next(),f.next(),f.next(),f.next(),f.next()]
        if i < N:        
            sample.append(read_info)
        elif random.random() < N/float(i+1):
            replace = random.randint(0,len(sample)-1)
            sample[replace] = read_info
        i+=1
    f.close()

f = open(sys.argv[2]+'/subsampled_shuffled.fastq','w')
for read_info in sample:
    f.write("".join(read_info))
f.close()