#!/usr/bin/python
import sys
import random
import glob

N = int(sys.argv[1])     #number of left reads required per sample for 10x total coverage; (10*genome size) / (read size*2*num_samples)
sample = [];
i=0

paired1=glob.glob(sys.argv[2]+"/*R1*fastq")
paired2=glob.glob(sys.argv[2]+"/*R2*fastq")
allfq=glob.glob(sys.argv[2]+"/*fastq")

for fq in paired1:
    print 'Subsampling '+fq
    f1=open(fq, 'r')
    f2=open(fq.replace('R1','R2'), 'r')
    
    for line in f1:
        if not line: break      #have less than desired coverage
        read_info=[line,f1.next(),f1.next(),f1.next(),f2.next(),f2.next(),f2.next(),f2.next()]
        if i < N:        
            sample.append(read_info)
        elif random.random() < N/float(i+1):
            replace = random.randint(0,len(sample)-1)
            sample[replace] = read_info
        i+=1
    f1.close()
    f2.close()
print(str(len(sample))+' pairs of reads sampled')

if not paired1:
    N=N*2
    for fq in allfq:
        if (fq not in paired1 and fq not in paired2):
            print 'Subsampling '+fq
            f1=open(fq, 'r')
            for line in f1:
                if not line: break      #have less than desired coverage
                read_info=[line,f1.next(),f1.next(),f1.next()]
                if i < N:        
                    sample.append(read_info)
                elif random.random() < N/float(i+1):
                    replace = random.randint(0,len(sample)-1)
                    sample[replace] = read_info
                i+=1
            f1.close()
            
    f = open(sys.argv[2]+'/subsampled.fastq','w')
    for read_info in sample:
        f.write("".join(read_info))
    f.close()
    
else:
    f = open(sys.argv[2]+'/subsampled_shuffled.fastq','w')
    for read_info in sample:
        f.write("".join(read_info))
    f.close()