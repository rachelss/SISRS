#!/usr/bin/env python3
#This script assumes a single pileup file per taxon
#Add an assert statement to test this
import os
import sys
from collections import Counter
from Bio import SeqIO
import glob

#get combined pileup info
def getallbases(path):
    assert len(glob.glob1(path,"*.pileups"))==1,'More than one pileup file in'+path
    allbases=dict()
    with open (path+'/'+os.path.basename(path)+'.pileups',"r") as filein:
        for line in filein:
            splitline=line.split()
            if len(splitline)>4:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                cleanBases=getCleanList(ref,bases)
                assert len(cleanBases) == int(num), 'bases are being counted incorrectly: '+ str(bases) + ' should have '+str(num)+' bases, but it is being converted to '+"".join(cleanBases)
                finalBase=getFinalBase_Specific(cleanBases)
                allbases[loc]=finalBase
    return allbases

def getCleanList(ref,bases):
    bases=bases.replace('.',ref) #insert ref base
    bases=bases.replace(',',ref)
    bases=bases.upper() #everything in uppercase
    bases=list(bases)

    okbases=['A','C','G','T','*']
    indels=['+','-']

    new_base_list=[]
    ibase_list = iter(bases)
    for b in ibase_list:
        if b in okbases:
            new_base_list.append(b)         #Get base
        elif b in indels:                   #skip indels
            i = int(next(ibase_list))
            j = str(next(ibase_list))
            if str.isdigit(j):
                skip=int(str(i)+j)
                while skip>0:
                    z=next(ibase_list)
                    skip=skip-1
            else:
                while i>1:
                    z=next(ibase_list)
                    i = i-1
        elif b=='^':                        #skip read qual noted at end of read
            z=next(ibase_list)

    return new_base_list

def getFinalBase_Specific(cleanBases):
    most_common = sorted(Counter(cleanBases).most_common())
    finalBase=(most_common[0][0])
    if finalBase == '*':
        finalBase = 'N'
    return finalBase

###############################################
#if __name__ == "__main__":
#    path=sys.argv[1]
#    contig_file = sys.argv[2]

def main(path, contig_file):

    allbases=getallbases(path)      #dictionary of combined pileups - locus/pos:bases(as list)
    if len(allbases)==0:
        print('No data for '+path)
        sys.exit(1)

    # Read contig fasta file into dictionary with sequence ID as the key
    contig_handle = open(contig_file, "r")
    fasta_seq = SeqIO.parse(contig_handle, 'fasta')
    fasta_dict = {read.id:list(str(read.seq)) for read in fasta_seq}
    contig_handle.close()

    for locus_pos,base in sorted(allbases.items()):
        locus,pos = locus_pos.split('/')
        fasta_dict[locus][int(pos)-1] = base

    #output = open(path+'/contigs.fa', 'wb')
    output = open(path+'/contigs.fa', 'w')
    for l,seq in sorted(fasta_dict.items()):
        output.write('>'+str(l)+"\n"+"".join(seq)+"\n")
    output.close()


if __name__ == '__main__':
    path = sys.argv[1]
    contig_file = sys.argv[2]
    main(path, contig_file)
