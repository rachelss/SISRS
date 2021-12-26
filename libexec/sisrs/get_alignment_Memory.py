#!/usr/bin/env python2
import os
import re
import glob
from collections import Counter,defaultdict
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from itertools import izip

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
        print str(len(self.locations))+' total variable sites (alignment.nex)'
        for i in range(len(self.locations)):
            bases = [self.species_data[sp][i] for sp in self.species_data if self.species_data[sp][i] in ['A','C','G','T','-']]     #bases for that site
            c = Counter(bases).most_common(5)
            if c[1][1]==1:
                self.single.append(1)
            else:
                self.single.append(0)
            self.flag.append(len(c))

        print str(self.single.count(1))+' variable sites are singletons'

        return self.flag.count(2)       # number of biallelic sites

def makerefdict(reffasta):
    ''' takes fasta file, returns dict of contig_name:[sequence]    '''

    filein=open(reffasta, "rU") #read fasta file
    refdict = SeqIO.to_dict(SeqIO.parse(filein, "fasta"))
    filein.close()

    for k in refdict:
        refdict[k]=list(str(refdict[k].seq))

    return refdict      #name:sequence_as_list

def sep_cigar(cigar):
    sepcigar = re.findall('[A-Z]+|[0-9]+',cigar)
    return sepcigar[::2], sepcigar[1::2]       #nums, code

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
                chro=nodedict[node].chr
                contig_map=nodedict[node].start
                flag = nodedict[node].flag
                cigar = nodedict[node].cigar
                if flag == '1':
                    pos=nodedict[node].length-pos+1
                pos = adjust_mapping(pos,cigar,0)
                pos=contig_map+pos-1   #1 based position
                alignment.ref_loc.append(chro+'_'+str(pos)) #add chr_pos
                if (pos-1) < len(ref[chro]):
                    refbase=ref[chro][(pos-1)] #chr in dictionary, site within seq-as-list   0 based position
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

def get_ref_info(alignment,reference,mainfolder,assembler):
    ''' #identify the chromosome and starting position for each contig '''

    nodedict = dict((loc.split('/')[0],'X') for loc in set(alignment.locations))    #make each contig name a key with an empty value (if it has a snp)
    if reference is 'X':
        ref=0
    else:
        ref=makerefdict(reference)        #make dict for reference chromosomes
        samfile=open(mainfolder+'/'+assembler+'output/align_contigs.sam','r') #sam file contains alignment of contigs to reference ref genome - has chr and pos
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
    return nodedict,ref

def write_alignment(fi,alignment,numbi):
    spp = sorted(alignment.species_data.keys())

    if len(alignment.ref) == 0:
        ntax = str(len(alignment.species_data))
    else:
        ntax = str(len(alignment.species_data)+1)

    ALIGNMENT=open(fi,'w')
    ALIGNMENTBI=open(fi.replace('.nex','_bi.nex'),'w')
    ALIGNMENTPI=open(fi.replace('.nex','_pi.nex'),'w')

    ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(alignment.locations))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENT.write('[ '+ " ".join(alignment.ref_loc)+' ]'+"\n")
    ALIGNMENT.write('[ '+ " ".join(alignment.locations)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENT.write(species+"\t"+"".join(alignment.species_data[species])+"\n")


    bi_ref_loc = [alignment.ref_loc[i] for i in range(len(alignment.locations)) if alignment.flag[i] == 2 and alignment.single[i] == 0]
    bi_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.flag[i] == 2 and alignment.single[i] == 0]
    bi_sp_data={}
    for species in spp:
        bi_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.flag[i] == 2 and alignment.single[i] == 0]
    if len(alignment.ref) > 0:
        bi_sp_data['reference'] = [alignment.ref[i] for i in range(len(alignment.locations)) if alignment.flag[i] == 2 and alignment.single[i] == 0]

    ALIGNMENTBI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(bi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTBI.write('[ '+ " ".join(bi_ref_loc)+' ]'+"\n")
    ALIGNMENTBI.write('[ '+ " ".join(bi_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTBI.write(species+"\t"+("".join(bi_sp_data[species]))+"\n")
    print str(len(bi_loc))+' total biallelic sites excluding singletons (alignment_bi.nex)'

    pi_ref_loc = [alignment.ref_loc[i] for i in range(len(alignment.locations)) if alignment.single[i] == 0]
    pi_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.single[i] == 0]
    pi_sp_data={}
    for species in spp:
        pi_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.single[i] == 0]
    if len(alignment.ref) > 0:
        pi_sp_data['reference'] = [alignment.ref[i] for i in range(len(alignment.locations)) if alignment.single[i] == 0]

    ALIGNMENTPI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(pi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTPI.write('[ '+ " ".join(pi_ref_loc)+' ]'+"\n")
    ALIGNMENTPI.write('[ '+ " ".join(pi_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTPI.write(species+"\t"+("".join(pi_sp_data[species]))+"\n")
    print str(len(pi_loc))+' total variable sites excluding singletons (alignment_pi.nex)'

    if len(alignment.ref) > 0:
        ALIGNMENT.write('reference'+"\t"+("".join(alignment.ref))+"\n")

    ALIGNMENT.write(';\nend;')
    ALIGNMENTBI.write(';\nend;')
    ALIGNMENTPI.write(';\nend;')

    ALIGNMENT.close()
    ALIGNMENTBI.close()
    ALIGNMENTPI.close()

def get_phy_sites(mainfolder,assembler,num_missing):

    #Fetch contig data
    contigList=glob.glob(mainfolder+'/' + assembler +'output/contigs_LocList')
    assert len(contigList) > 0, 'Total site list not found in assembly folder'

    #Fetch sorted species data
    dataLists = sorted(glob.glob(mainfolder+'/*/*_LocList'))
    dataLists.remove(mainfolder+'/' + assembler +'output/contigs_LocList')
    splist=[os.path.basename(os.path.dirname(path)) for path in dataLists]
    speciesCount=len(dataLists)
    assert len(dataLists) > 0, 'No species had data from the pileup'

    allLists = contigList+dataLists

    alignment = Alignment()
    alignment.species_data = {species: [] for species in splist}

    files = [open(i, "r") for i in allLists]
    for rows in izip(*files):
        rowList = map(lambda foo: foo.replace('\n', ''), list(rows))
        speciesData = rowList[1:(speciesCount+1)]
        if speciesData.count("N")<=num_missing and len(set(filter(lambda a: a != "N", speciesData)))>1:
            alignment.locations.append(rowList[0])
            for j in range(0,(speciesCount)):
                alignment.species_data[splist[j]].append(speciesData[j])
    return alignment

num_missing=int(sys.argv[1])
reference = sys.argv[2]
mainfolder = sys.argv[3]
assembler = sys.argv[4]
basecomplement = {'a':'t', 'c':'g', 't':'a', 'g':'c', 'A':'t', 'C':'g', 'T':'a', 'G':'c','-':'-'}

alignment=get_phy_sites(mainfolder,assembler,num_missing)
numbi=alignment.numsnps()
nodedict,ref = get_ref_info(alignment,reference,mainfolder,assembler)
alignment = add_ref_info(alignment,ref,nodedict)
alignment = write_alignment(mainfolder+'/alignment.nex',alignment,numbi)
