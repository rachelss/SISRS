#!/usr/bin/env python2
import os
import sys
##############

def createPosList(path,assembler):
    posList=[]
    with open(path+"/"+assembler+"output/contigs_SeqLength.tsv","r") as filein:
        for line in iter(filein):
            splitline=line.split()
            lengthList=range(1,(int(splitline[1])+1))
            for x in lengthList:
                posList.append(splitline[0] +'/'+str(x))
    printList = open(path+'/'+assembler+'output/contigs_LocList','w')
    for item in posList:
        print>>printList,item
    printList.close()
    sys.stdout.write("Site list created: " + str(len(posList)) + " total sites\n")

######################################################

path = sys.argv[1]
assembler = sys.argv[2]
createPosList(path,assembler)
