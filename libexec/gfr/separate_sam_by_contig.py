#! /usr/local/bin/python
import sys

folder=sys.argv[1].split('/')[0]
samfile=open(sys.argv[1],'r')   #open file
nodedict=dict()
for line in samfile:      #go through file
    if line.startswith('@'):    #skip header
        continue
    else:
        splitline=line.split()
        node = splitline[2]     #get contig aligned to
        if node is not '*':
            if node not in nodedict:
                nodedict[node]=list()
            nodedict[node].append(line)
samfile.close()

for k,v in nodedict.iteritems():
    outfile = open(folder+'/'+k+'.sam','a')
    for item in v:
        outfile.write("%s" % item)
    outfile.close()