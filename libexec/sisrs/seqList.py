#!/usr/bin/env python2
import os
import sys
import pandas as pd
##############

def createPosList(path,assembler):
    #Read in sequence lengths from SISRS output in <PATH>/<assembler>output/contigs_SeqLength.tsv; Sort by contig ID
    df = pd.read_csv(path+"/"+assembler+"output/contigs_SeqLength.tsv",sep="\t",header=None)
    df.sort_values([0],inplace=True)
    df=df.reset_index(drop=True)

    #Fill in list
    j=0
    posList = []
    for row in df.iterrows():
        contig=row[1][0]
        length=row[1][1]
        for i in range(1,length+1):
            pos = contig + "/" + str(i)
            posList.insert(j,pos)
            j+=1

    #Write file to <PATH>/<assembler>output/contigs_PosList
    printList = open(path+"/"+assembler+"output/contigs_PosList", 'w')
    for item in posList:
        print>>printList, item
    printList.close()
    print "Site list created: " + str(len(posList)) + " total sites"
    sys.stdout.flush()


if __name__ == "__main__":

    #Read in arguments
    path = sys.argv[1]
    assembler = sys.argv[2]
    createPosList(path,assembler)
