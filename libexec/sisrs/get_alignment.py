#!/usr/bin/env python2
import os
import re
import cPickle
from collections import Counter
import sys

#########################
def makerefdict(reffasta):
    refdict=dict()
    filein=open(reffasta, 'r') #read fasta file
    for line in filein.xreadlines(): #go through each line
        if line.startswith('>'): #get chr name
            splitline=line.split()
            refdict[splitline[0].replace('>', '')]=[]
        else:
            refdict[splitline[0].replace('>', '')].append(line.rstrip()) #get seq
    filein.close()
    for k in refdict.iterkeys():
        refdict[k]="".join(refdict[k])
        refdict[k]=list(refdict[k])
#        print 'chr: '+k+"\t"+'length: '+str(len(refdict[k]))
    return refdict #name:sequence as list
#########################
#make dict for reference chromosomes
if sys.argv[2] is not 'X':
    ref=makerefdict(sys.argv[2])

pathlist=[]
allbases=dict()

#list of directories that have pileup files
for path,dirs,files in os.walk(sys.argv[3]):
    for file in files:
        if file.endswith(".pkl"):
            pathlist.append(path)
pathlist=list(set(pathlist))
pathlist.sort()
    
alllocs=[]
for species in pathlist:
    # read python dict back from the file
    print 'Reading data: '+species
    pkl_file = open(species+'/pruned_dict.pkl', 'rb')
    sp_bases = cPickle.load(pkl_file)
    pkl_file.close()
    for key in sp_bases.keys():
        alllocs.append(key)
    allbases[species]=sp_bases

#get sorted list of unique position names
alllocs=list(set(alllocs))
alllocs.sort()  
    
aligndic=dict()
aligndic['locations']=[] #list of locations - keep in same order as list of snps
for species in pathlist:
    aligndic[species]=[]
    
aligndic2=dict()
aligndic2['locations']=[] #list of locations - keep in same order as list of snps - biallelic data only
for species in pathlist:
    aligndic2[species]=[]

numsnps=0
numbi=0
for loc in alllocs: #go through each location
    flag=0
    snp=[]
    for species in pathlist: #get list of bases for this location for all species
        if loc in allbases[species]:
#            print loc,allbases[species][loc]
            snp.append(allbases[species][loc])
    snp=list(("".join(snp)).replace('N','')) #make list into string, replace N's, put back into list
    basecounts=Counter(snp).most_common()
    if len(snp)<(len(pathlist)-int(sys.argv[1])): #missing info for more than X indivs
        flag=2
#        print snp
    elif(len(basecounts))<2: #no variation or empty - flag
        flag=2
#    elif basecounts[1][1]==1: #if there are two bases and one's a singleton OR three bases and two are singletons etc - flag
#        flag=2
    elif(len(basecounts))>2: #not biallelic - flag for full dict
        flag=1
    if flag<2: #if there is variation and it's not due to a singleton
        print loc+' is a SNP'
        numsnps+=1
        aligndic['locations'].append(loc) #add position name to list
        for species in pathlist: #add base to the list for that species (or add 'N')
            if loc in allbases[species]:
                aligndic[species].append(allbases[species][loc])
            else:
                aligndic[species].append('N')
        if flag==0:
            print 'biallelic'
            numbi+=1
            aligndic2['locations'].append(loc) #add position name to list
            for species in pathlist: #add base to the list for that species (or add 'N')
                if loc in allbases[species]:
                    aligndic2[species].append(allbases[species][loc])
                else:
                    aligndic2[species].append('N')
        print snp
print str(numsnps)+' snps'
print str(numbi)+' biallelic snps'

#identify the chromosome and starting position for each contig
nodedict=dict()
nodes=list(set(aligndic['locations'])) #list of contigs
for node in nodes:
    node,site=node.split('/')
    nodedict[node]='X' #make each contig name a key with an empty value

if sys.argv[2] is not 'X':
    samfile=open(sys.argv[3]+'/velvetoutput/align_contigs.sam','r') #sam file contains alignment of contigs to reference ref genome - has chr and pos
    for line in samfile.xreadlines():
        if line.startswith('NODE'):
            splitline=line.split()
            node=splitline[0] #contig name
            chr=splitline[2]
            site=int(splitline[3])
            if node in nodedict:
                if site>0:
                    nodedict[node]=chr,site #assign tuple containing chr,site as value to key contigname

aligndic['reference']=[]
aligndic2['reference']=[]

ALIGNMENT=open(sys.argv[3]+'/alignment.nex','w')
ALIGNMENT2=open(sys.argv[3]+'/alignment_bi.nex','w')

ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+str((len(pathlist))+1)+' NCHAR='+str(len(aligndic['locations']))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n[')
for loc in aligndic['locations']: #write list of locations in order - translate into chr,site info
    node,site=loc.split('/')
    if nodedict[node] is not 'X': #contig maps to reference
        chr=nodedict[node][0]
#        print chr
        pos=nodedict[node][1]+int(site)-1   #1 based position
#        print pos
        loc2=chr+'_'+str(pos) #chr_pos
        refbase=ref[chr][(pos-1)] #chr in dictionary, site within seq-as-list   0 based position
        aligndic['reference'].append(refbase.upper())
    else:
        loc2='X'
        aligndic['reference'].append('N') #if we don't know where the contig is mapped, we use N for the ref base
    ALIGNMENT.write(loc2+' ')
ALIGNMENT.write(']\n')
for species in pathlist: #write sequences for each species
    ALIGNMENT.write(species+"\t"+("".join(aligndic[species]))+"\n")
ALIGNMENT.write('reference'+"\t"+("".join(aligndic['reference']))+"\n")
ALIGNMENT.write(';\nend;')
ALIGNMENT.close()

ALIGNMENT2.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+str((len(pathlist))+1)+' NCHAR='+str(len(aligndic2['locations']))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n[')
for loc in aligndic2['locations']: #write list of locations in order - translate into chr,site info
    node,site=loc.split('/')
    if nodedict[node] is not 'X': #contig maps to reference
        chr=nodedict[node][0]
        pos=nodedict[node][1]+int(site)-1
        loc2=chr+'_'+str(pos) #chr_pos
        refbase=ref[chr][(pos-1)] #chr in dictionary, site within seq-as-list
        aligndic2['reference'].append(refbase.upper())
    else:
        loc2='X'
        aligndic2['reference'].append('N') #if we don't know where the contig is mapped, we use N for the ref base
    ALIGNMENT2.write(loc2+' ')
ALIGNMENT2.write(']\n')
for species in pathlist: #write sequences for each species
    ALIGNMENT2.write(species+"\t"+("".join(aligndic2[species]))+"\n")
ALIGNMENT2.write('reference'+"\t"+("".join(aligndic2['reference']))+"\n")
ALIGNMENT2.write(';\nend;')
ALIGNMENT2.close()