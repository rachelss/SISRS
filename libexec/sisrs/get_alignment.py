#!/usr/bin/env python2
import os
import re
import cPickle
from collections import Counter
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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
        
class Alignment:
    def __init__(self,locations=[],species_data=dict(),ref=[],ref_loc=[],flag=[],single=[]):
        self.locations = locations
        self.species_data = species_data
        self.ref = ref
        self.ref_loc = ref_loc
        self.flag = flag
        self.single = single
        
    def numsnps(self):
        print str(len(self.locations))+' variable sites'
        singletons,bi=0,0
        for i in range(len(self.locations)):
            bases = [self.species_data[sp][i] for sp in self.species_data if self.species_data[sp][i] in ['A','C','G','T']]
            c = Counter(bases).most_common(4)
            if c[1][1]==1:
                singletons+=1
                self.single.append(1)
            else:
                self.single.append(0)
            if len(c) == 2:
                bi+=1
            self.flag.append(len(c))
                
        print str(bi)+' biallelic sites'
        print str(singletons)+' singletons'
        
        return bi

def makerefdict(reffasta):
    filein=open(reffasta, "rU") #read fasta file
    refdict = SeqIO.to_dict(SeqIO.parse(filein, "fasta"))
    filein.close()

    for k in refdict:
        refdict[k]=list(str(refdict[k].seq))
        
    return refdict #name:sequence_as_list

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
        for key in sp_bases:
            alllocs.append(key)
        allbases[species]=sp_bases
    
    alllocs=list(set(alllocs))
    alllocs.sort()  
        
    return allbases,alllocs

def get_phy_sites(pathlist,allbases,alllocs,num_missing):
    alignment = Alignment()
    for species in pathlist:
        alignment.species_data[species]=[]
    for loc in alllocs: #go through each location
        snp = [allbases[species][loc] for species in pathlist if loc in allbases[species]]
        c = Counter(snp)
        if (len(pathlist)-len(snp))+c['N'] > num_missing:     #too many missing, go to next loc
            continue
        del c['N']
        if len(c) == 1:     #no variation
            continue
        else:
            alignment.locations.append(loc)
            for species in pathlist: #add base to the list for that species (or add 'N')
                if loc in allbases[species]:
                    alignment.species_data[species].append(allbases[species][loc])
                else:
                    alignment.species_data[species].append('N')

    return alignment

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

def add_ref_info(alignment,ref,nodedict):
    if ref == 0:
        alignment.ref_loc = ['X'] * len(alignment.locations)        #no reference, just use X for all ref locations
    else:
        for loc in alignment.locations:
            node,site=loc.split('/')
            pos=int(site)
            if nodedict[node] is not 'X': #contig maps to reference
                chr=nodedict[node].chr
                contig_map=nodedict[node].start
                flag = nodedict[node].flag
                cigar = nodedict[node].cigar
                if flag == '1':
                    pos=nodedict[node].length-pos+1
                pos = adjust_mapping(pos,cigar,0)
                pos=contig_map+pos-1   #1 based position
                alignment.ref_loc.append(chr+'_'+str(pos)) #add chr_pos
                if (pos-1) < len(ref[chr]):
                    refbase=ref[chr][(pos-1)] #chr in dictionary, site within seq-as-list   0 based position
                    if flag == '1':
                        if refbase in basecomplement:
                            refbase = basecomplement[refbase]
                        else:
                            refbase = 'N'
                    alignment.ref.append(refbase.upper()) #add ref base
            else:
                alignment.ref_loc.append('X')
                alignment.ref.append('N')

    return alignment

def write_alignment(fi,alignment,numbi):
    if len(alignment.ref) == 0:
        ntax = str(len(alignment.species_data))
    else:
        ntax = str(len(alignment.species_data)+1)
    
    ALIGNMENT=open(fi,'w')
    ALIGNMENTBI=open(fi.replace('.nex','_bi.nex'),'w')
    ALIGNMENTPI=open(fi.replace('.nex','_pi.nex'),'w')
    ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(alignment.locations))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTBI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(numbi)+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
        
    ALIGNMENT.write('[ '+ " ".join(alignment.ref_loc)+' ]'+"\n")
    ALIGNMENT.write('[ '+ " ".join(alignment.locations)+' ]'+"\n")
    for species in pathlist: #write sequences for each species
        ALIGNMENT.write(species.replace('./','')+"\t"+("".join(alignment.species_data[species]))+"\n")
        
    bi_ref_loc,bi_loc,bi_ref=[],[],[]
    pi_ref_loc,pi_loc,pi_ref=[],[],[]
    bi_sp_data,pi_sp_data = dict(),dict()
    for species in pathlist:
        bi_sp_data[species] = []
        pi_sp_data[species] = []
    for i in range(len(alignment.locations)):
        if alignment.flag[i] == 2:
            bi_ref_loc.append(alignment.ref_loc[i])
            bi_loc.append(alignment.locations[i])
            for species in pathlist:
                bi_sp_data[species].append(alignment.species_data[species][i])
            if len(alignment.ref) > 0:
                bi_ref.append(alignment.ref[i])
        if alignment.single[i] == 0:
            pi_ref_loc.append(alignment.ref_loc[i])
            pi_loc.append(alignment.locations[i])
            for species in pathlist:
                pi_sp_data[species].append(alignment.species_data[species][i])
            if len(alignment.ref) > 0:
                pi_ref.append(alignment.ref[i])
                
    ALIGNMENTBI.write('[ '+ " ".join(bi_ref_loc)+' ]'+"\n")
    ALIGNMENTBI.write('[ '+ " ".join(bi_loc)+' ]'+"\n")
    for species in pathlist: #write sequences for each species
        ALIGNMENTBI.write(species.replace('./','')+"\t"+("".join(bi_sp_data[species]))+"\n")
    
    ALIGNMENTPI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(pi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')    
    ALIGNMENTPI.write('[ '+ " ".join(pi_ref_loc)+' ]'+"\n")
    ALIGNMENTPI.write('[ '+ " ".join(pi_loc)+' ]'+"\n")
    for species in pathlist: #write sequences for each species
        ALIGNMENTPI.write(species.replace('./','')+"\t"+("".join(pi_sp_data[species]))+"\n")
        
    if len(alignment.ref) > 0:
        ALIGNMENT.write('reference'+"\t"+("".join(alignment.ref))+"\n")
        ALIGNMENTBI.write('reference'+"\t"+("".join(bi_ref))+"\n")
        ALIGNMENTPI.write('reference'+"\t"+("".join(pi_ref))+"\n")
    ALIGNMENT.write(';\nend;')
    ALIGNMENTBI.write(';\nend;')
    ALIGNMENTPI.write(';\nend;')
    ALIGNMENT.close()
    ALIGNMENTBI.close()
    ALIGNMENTPI.close()
    
    return 1

#########################
num_missing = int(sys.argv[1])
basecomplement = {'a':'t', 'c':'g', 't':'a', 'g':'c', 'A':'t', 'C':'g', 'T':'a', 'G':'c'}

pathlist = get_pileup_files(sys.argv[3])
num_species = len(pathlist)
allbases,alllocs = read_pkls(pathlist)      #dict of species:(loc:base), sorted list of unique position names

alignment=get_phy_sites(pathlist,allbases,alllocs,num_missing)       #return Alignment object of informative sites
numbi = alignment.numsnps()     #prints numbers of snps, biallelic snps, and singletons

#identify the chromosome and starting position for each contig
nodedict = dict((loc.split('/')[0],'X') for loc in set(alignment.locations))    #make each contig name a key with an empty value (if it has a snp)
if sys.argv[2] is not 'X':
    ref=makerefdict(sys.argv[2])        #make dict for reference chromosomes
    samfile=open(sys.argv[3]+'/velvetoutput/align_contigs.sam','r') #sam file contains alignment of contigs to reference ref genome - has chr and pos
    for line in samfile:
        if not line.startswith('@'):    #skip header
            splitline=line.split()
            if splitline[0] in nodedict:        #only interested if contig has a snp
                flag = bin(int(splitline[1]))
                flagl = flag.split('b')
                if len(flagl[1])>=3:
                    if flagl[1][-3] == '1':
                        continue                #if contig is unmapped, continue
                if len(flagl[1])>=5:
                    flag = flagl[1][-5]     #check for reverse mapping
                else:
                    flag = '0'    

                chro=splitline[2]
                start=int(splitline[3])
                cigar=(splitline[5])
                end = start + (adjust_mapping(len(splitline[9]),cigar,0)) - 1
                nodedict[splitline[0]]=Scaffold(chro,start,end,len(splitline[9]),flag,cigar)
else:
    ref=0

alignment = add_ref_info(alignment,ref,nodedict)
alignment = write_alignment(sys.argv[3]+'/alignment.nex',alignment,numbi)