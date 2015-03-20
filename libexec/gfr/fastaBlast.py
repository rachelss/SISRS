#! /usr/local/bin/python
#modified from Bienvenido Velez
import sys
import os.path
from Bio import SeqIO,GenBank,Entrez
from Bio.Blast import NCBIXML,NCBIWWW

def runBlastn(fastaSequence,outputFileName):
    'Conducts a nucleotide Blast run on the given fastaSequence and save the blast run output in XML format'
    result_handle = NCBIWWW.qblast("blastn", "nr",fastaSequence)
    # save the blast results for later, in case we want to look at them
    save_file = open(outputFileName, "w")
    blast_results = result_handle.read()
    save_file.write(blast_results)
    save_file.close()
    return blast_results    

###########################
Entrez.email = 'Rachel.Schwartz@asu.edu'

all_contigs_file = open("velvetoutput/contigs.fa", "rU")
allcontigs = SeqIO.to_dict(SeqIO.parse(all_contigs_file, "fasta"))
all_contigs_file.close()

locus_list=open('loci.txt','r')
loci=locus_list.readlines()
locus_list.close()

for l in loci:
    if not os.path.exists('loci/'+l.rstrip()+'_extra.fa'):
        blast_results=runBlastn(allcontigs[l.rstrip()].seq,l+'.out')
        blast_out = open(l+'.out')
        blast_record = NCBIXML.read(blast_out)
        
        a={}
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.1:
                    a[alignment.hit_id]=hsp.sbjct
        
        fileout=open('loci/'+l.rstrip()+'_extra.fa','w')
        sp=[]
        for k in a.keys():
            handle = Entrez.efetch(db="nucleotide", id=k, retmode="xml")
            records = Entrez.read(handle)
            if '_'.join(records[0]['GBSeq_source'].split()[0:2]) not in sp:
                fileout.write('>'+'_'.join(records[0]['GBSeq_source'].split()[0:2])+"\n")
                fileout.write(a[k]+"\n")
                sp.append('_'.join(records[0]['GBSeq_source'].split()[0:2]))
    
        fileout.close()

#while read i; do
#    F=$( echo $i | sed 's/\.[^.]*$//' ).fa
#    cat $F ${i}_extra.fa > ${F}_all.fa
#done < loci.txt
#find loci/ -type f -name '*_all.fa' | parallel -j "${PROCESSORS}" "mafft {} > {.}_align.fa"