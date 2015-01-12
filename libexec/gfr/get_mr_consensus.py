#! /usr/local/bin/python
#take sam file aligned to single contig and produce whole alignment
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter
import glob

#######################
def get_consensus(data):
    pos_bases_counts = Counter(data).most_common()
    c_base = pos_bases_counts[0][0]
    return c_base

######################
bases = ['A','C','G','T','a','c','g','t']

genes = glob.glob("*_aligned.fasta")
geneseqs = []

for gene in genes:
    final_seq=[]
    handle = open(gene, "rU")
    allseqs = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    
    for i in range(len(allseqs[0])):
        data = [seq[i] for seq in allseqs if seq[i] in bases]
        if len(data)>0:
            b = get_consensus(data)
            final_seq.append(b)
    
    final_seq_r = SeqRecord(Seq(''.join(final_seq)), id=gene.replace('_aligned.fasta',''))
    geneseqs.append(final_seq_r)

SeqIO.write(geneseqs,'ref_genes.fa', "fasta")