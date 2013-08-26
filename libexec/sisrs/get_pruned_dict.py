#!/usr/bin/env python2
import os
import sys
import cPickle

allbases=dict()
path=sys.argv[1]
minread=sys.argv[2]
minread=int(minread)
for fi in os.listdir(path):
    if fi.endswith("pileups"):
#        print fi
        filein=open(path+'/'+fi,'r')
        for line in filein:
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                bases=bases.replace('.',ref) #insert ref base
                bases=bases.replace(',',ref)
                bases=bases.upper() #everything in uppercase
                bases=list(bases)
                loc=node+'/'+pos
                if loc in allbases:
                    allbases[loc].append(bases)
                else:
                    allbases[loc]=bases
        filein.close()

#prune by whether there's enough info and species are fixed
for loc in allbases.iterkeys():
    if((allbases[loc].count('A')+allbases[loc].count('C')+allbases[loc].count('G')+allbases[loc].count('T'))< minread): #not enough info
        allbases[loc]='N'
    else:
        if(bool(allbases[loc].count('A'))+bool(allbases[loc].count('C'))+bool(allbases[loc].count('G'))+bool(allbases[loc].count('T')))>1: #not fixed
            allbases[loc]='N'
        else: #set base for this species
            if(allbases[loc].count('A')>0):
                allbases[loc]='A'
            elif(allbases[loc].count('C')>0):
                allbases[loc]='C'
            elif(allbases[loc].count('G')>0):
                allbases[loc]='G'
            elif(allbases[loc].count('T')>0):
                allbases[loc]='T'
#    print allbases[loc]

output = open(path+'/pruned_dict.pkl', 'wb')
cPickle.dump(allbases, output, cPickle.HIGHEST_PROTOCOL)
output.close()