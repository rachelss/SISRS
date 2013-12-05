#!/usr/bin/env python2
import os
import re
import cPickle
from collections import Counter
import sys

#########################
class Scaffold:
    def __init__(self,chr,start,end,length=None,flag=None,cigar=None):
        self.chr = chr
        self.start = start
        self.end = end
        self.length = length
        self.flag = flag
        self.cigar = cigar
        
class Loc:
    def __init__(self,scaff_loc,flag):
        self.scaff_loc = scaff_loc
        self.flag = flag

def makerefdict(reffasta):
    refdict=dict()
    filein=open(reffasta, 'r') #read fasta file
    for line in filein: #go through each line
        if line.startswith('>'): #get chr name
            splitline=line.split()
            refdict[splitline[0].replace('>', '')]=[]
        else:
            refdict[splitline[0].replace('>', '')].append(line.rstrip()) #get seq
    filein.close()
    for k in refdict.iterkeys():
        refdict[k]="".join(refdict[k])
        refdict[k]=list(refdict[k])
    return refdict #name:sequence as list

def get_pileup_files(f):
    pathlist=[]
    for path,dirs,files in os.walk(f):
        for file in files:
            if file.endswith(".pkl"):
                pathlist.append(path)
    pathlist=list(set(pathlist))
    pathlist.sort()
    
    return pathlist

def read_pkls(pathlist):
    allbases=dict()
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
    
    alllocs=list(set(alllocs))
    alllocs.sort()  
        
    return allbases,alllocs

def TestSNP(loc):
    flag=0
    snp=[]
    for species in pathlist: #get list of bases for this location for all species
        if loc in allbases[species]:
            snp.append(allbases[species][loc])
    snp=list(("".join(snp)).replace('N','')) #make list into string, replace N's, put back into list
    basecounts=Counter(snp).most_common()
    if len(snp)<(num_species-num_missing): #missing info for more than X indivs
        flag=2
    elif(len(basecounts))<2: #no variation or empty - flag
        flag=2
    elif(len(basecounts))>2: #not biallelic - flag for full dict
        flag=1
        
    return flag

def make_aligndic(pathlist):
    aligndic=dict()
    aligndic['locations']=[] #list of locations - keep in same order as list of snps
    aligndic['flags'] = []
    for species in pathlist:
        aligndic[species]=[]
    aligndic['reference']=[]
        
    return aligndic

def sep_cigar(cigar):
    d=list()
    code=list()
    nums=list()
    for i,c in enumerate(cigar):
        if c.isdigit():
            d.append(c)
        else:
            code.append(c)
            d="".join(d)
            nums.append(int(d))
            d=list()

    return nums, code

def adjust_mapping(readmap,cigar,flag):
    
    nums,code = sep_cigar(cigar)
    
    i=1
    newreadmap = readmap
    
    if flag == 1:
        code.reverse()
        nums.reverse()
    
    for j,num in enumerate(nums):

        if code[j] == 'S' or code[j] == 'I':
            newreadmap = newreadmap - num
        elif code[j] == 'D':
            newreadmap = newreadmap + num
            
        i = i+num
        if i>readmap:
            break
            
    return newreadmap

def write_alignment(fi,nchar):
    ALIGNMENT=open(fi,'w')
    ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(nchar)+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n[')
    for loc in aligndic['locations']: #write list of locations in order - translate into chr,site info
        node,site=loc.split('/')
        pos=int(site)
        if nodedict[node] is not 'X': #contig maps to reference
            chr=nodedict[node].chr
            contig_map=nodedict[node].start
    #        print chr
            flag = nodedict[node].flag
            cigar = nodedict[node].cigar
            if flag == '1':
                pos=nodedict[node].length-pos+1
            pos = adjust_mapping(pos,cigar,0)
            pos=contig_map+pos-1   #1 based position
    #           print pos
            loc2=chr+'_'+str(pos) #chr_pos
            if (pos-1) < len(ref[chr]):
                refbase=ref[chr][(pos-1)] #chr in dictionary, site within seq-as-list   0 based position
                if flag == '1':
                    if refbase in basecomplement:
                        refbase = basecomplement[refbase]
                    else:
                        refbase = 'N'
                aligndic['reference'].append(refbase.upper())
            else:
                loc2='X'
                aligndic['reference'].append('N')
        else:
            loc2='X'
            aligndic['reference'].append('N') #if we don't know where the contig is mapped, we use N for the ref base
        ALIGNMENT.write(loc2+' ')
    ALIGNMENT.write(']\n')
    ALIGNMENT.write('[ '+ " ".join(aligndic['locations'])+' ]'+"\n")
    for species in pathlist: #write sequences for each species
        ALIGNMENT.write(species.replace('./','')+"\t"+("".join(aligndic[species]))+"\n")
    if sys.argv[2] is not 'X':
        ALIGNMENT.write('reference'+"\t"+("".join(aligndic['reference']))+"\n")
    ALIGNMENT.write(';\nend;')
    ALIGNMENT.close()
    
    return 1
#########################
num_missing = int(sys.argv[1])
basecomplement = {'a':'t', 'c':'g', 't':'a', 'g':'c', 'A':'t', 'C':'g', 'T':'a', 'G':'c'}

#make dict for reference chromosomes
if sys.argv[2] is not 'X':
    ref=makerefdict(sys.argv[2])
else:
    ref=0

pathlist = get_pileup_files(sys.argv[3])
num_species = len(pathlist)
allbases,alllocs = read_pkls(pathlist)      #dict of species:(loc:base), sorted list of unique position names
aligndic = make_aligndic(pathlist)
aligndic2 = make_aligndic(pathlist)

numsnps,numbi=0,0
for loc in alllocs: #go through each location
    flag = TestSNP(loc)
    if flag<2: #if there is variation
        numsnps+=1
        aligndic['locations'].append(loc) #add position name to list
        aligndic['flags'].append(flag) #add position name to list
        for species in pathlist: #add base to the list for that species (or add 'N')
            if loc in allbases[species]:
                aligndic[species].append(allbases[species][loc])
            else:
                aligndic[species].append('N')
        if flag==0:
            numbi+=1
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
    for line in samfile:
        if line.startswith('@'):            #skip header
            continue
        else:
            splitline=line.split()            
            flag = bin(int(splitline[1]))
            flagl = flag.split('b')
            if len(flagl[1])>=3:
                if flagl[1][-3] == '1':
                    continue                #if scaffold is unmapped, continue
            if len(flagl[1])>=5:
                flag = flagl[1][-5]     #check for reverse mapping
            else:
                flag = '0'
                
            node=splitline[0] #contig name
            if node in nodedict:
                chr=splitline[2]
                start=int(splitline[3])
                cigar=(splitline[5])
                end = start + (adjust_mapping(len(splitline[9]),cigar,0)) - 1
                nodedict[node]=Scaffold(chr,start,end,len(splitline[9]),flag,cigar)
    ntax=str(len(pathlist)+1)
else:
    ntax=str(len(pathlist))

write_alignment(sys.argv[3]+'/alignment.nex',len(aligndic['locations']))
write_alignment(sys.argv[3]+'/alignment_bi.nex',sum(aligndic['flags']))