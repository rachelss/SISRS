#!/usr/bin/env python3

"""
    output alignment of sites useful for phylogenetics (nexus format)
    can specify max number of species missing for each site

    arguments:
        num_missing -- the number of species allowed to be missing at a site so that site will show in output alignment
        mainfolder -- folder that includes folders containing data
        assembler -- assembler used to get composite reference

    output:
        alignment.nex : nexus formatted alignment including position in composite reference genome; each site as up to num_missing missing data
        alignment_bi.nex : above but only biallelic sites
        alignment_pi.nex : above but only phylogenetically informative sites (no singletons)

"""

import os
import re
import glob
from collections import Counter,defaultdict
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#########################
class Loc:
    def __init__(self,scaff_loc,flag):
        self.scaff_loc = scaff_loc
        self.flag = flag

class Alignment:
    def __init__(self,locations=[],species_data=dict(),flag=[],single=[]):
        self.locations = locations
        self.species_data = species_data
        self.flag = flag
        self.single = single

    def numsnps(self):
        print(str(len(self.locations))+' total variable sites (alignment.nex)')
        for i in range(len(self.locations)):
            bases = [self.species_data[sp][i] for sp in self.species_data if self.species_data[sp][i] in ['A','C','G','T','-']]     #bases for that site
            c = Counter(bases).most_common(5)
            if c[1][1]==1:
                self.single.append(1)
            else:
                self.single.append(0)
            self.flag.append(len(c))

        print(str(self.single.count(1))+' variable sites are singletons')

        return self.flag.count(2)       # number of biallelic sites


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
    for rows in zip(*files):
        rowList = map(lambda foo: foo.replace('\n', ''), list(rows))
        speciesData = rowList[1:(speciesCount+1)]
        if speciesData.count("N")<=num_missing and len(set(filter(lambda a: a != "N", speciesData)))>1:
            alignment.locations.append(rowList[0])
            for j in range(0,(speciesCount)):
                alignment.species_data[splist[j]].append(speciesData[j])
    return alignment

def write_alignment(fi,alignment,numbi):
    spp = sorted(alignment.species_data.keys())
    ntax = str(len(alignment.species_data))

    #Process alignment.nex
    ALIGNMENT=open(fi,'w')
    ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(alignment.locations))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENT.write('[ '+ " ".join(alignment.locations)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENT.write(species+"\t"+"".join(alignment.species_data[species])+"\n")
    ALIGNMENT.write(';\nend;')
    ALIGNMENT.close()

    #Process alignment_bi.nex
    ALIGNMENTBI=open(fi.replace('.nex','_bi.nex'),'w')
    bi_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.flag[i] == 2 and alignment.single[i] == 0]
    bi_sp_data={}
    for species in spp:
        bi_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.flag[i] == 2 and alignment.single[i] == 0]

    ALIGNMENTBI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(bi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTBI.write('[ '+ " ".join(bi_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTBI.write(species+"\t"+("".join(bi_sp_data[species]))+"\n")
    print(str(len(bi_loc))+' total biallelic sites excluding singletons (alignment_bi.nex)')
    ALIGNMENTBI.write(';\nend;')
    ALIGNMENTBI.close()

    #Process alignment_pi.nex
    ALIGNMENTPI=open(fi.replace('.nex','_pi.nex'),'w')
    pi_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.single[i] == 0]
    pi_sp_data={}
    for species in spp:
        pi_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.single[i] == 0]

    ALIGNMENTPI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(pi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTPI.write('[ '+ " ".join(pi_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTPI.write(species+"\t"+("".join(pi_sp_data[species]))+"\n")
    print(str(len(pi_loc))+' total variable sites excluding singletons (alignment_pi.nex)')
    ALIGNMENTPI.write(';\nend;')
    ALIGNMENTPI.close()

#########################

def main(num_missing, mainfolder, assembler):

    alignment=get_phy_sites(mainfolder,assembler,num_missing)
    numbi=alignment.numsnps() #prints numbers of snps, biallelic snps, and singletons
    alignment = write_alignment(mainfolder+'/alignment.nex',alignment,numbi)

if __name__ == '__main__':
    num_missing = int(sys.argv[1])
    mainfolder = sys.argv[2]
    assembler = sys.argv[3]
    main(num_missing, mainfolder, assembler)
