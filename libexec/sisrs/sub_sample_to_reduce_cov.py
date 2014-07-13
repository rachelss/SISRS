#!/usr/bin/python
import sys
import random
import glob

N = int(sys.argv[1])     #final coverage
G = int(sys.argv[2])    #genome size
R = G*N/200    #number pairs of reads
sample = [];
i=0

for fq in glob.glob(sys.argv[3]+"/*shuffled*fastq"):
    print 'Subsampling '+fq
    f=open(fq, 'r')
    for line in f:
        if not line: break      #have less than desired coverage
        read_info=[line,f.next(),f.next(),f.next(),f.next(),f.next(),f.next(),f.next()]
        if i < R:        
            sample.append(read_info)
        elif random.random() < R/float(i+1):
            replace = random.randint(0,len(sample)-1)
            sample[replace] = read_info
        i+=1
    f.close()

out1 = open(sys.argv[3]+'/subsampled_R1.fastq','w')
out2 = open(sys.argv[3]+'/subsampled_R2.fastq','w')
for read_info in sample:
    out1.write("".join(read_info[:4]))
    out2.write("".join(read_info[4:]))
out1.close()
out2.close()