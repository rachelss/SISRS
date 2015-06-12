#!/usr/bin/env python2
from Bio import AlignIO
import sys
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

numsp=float(sys.argv[2])
bases=['a','c','g','t','A','C','G','T']
f=sys.argv[1]
fformat = sys.argv[1].split('.')[-1]
align = AlignIO.read(f,fformat)

for i in range(len(align[0])):      #go through each base in alignment
    info=[a for a in list(align[:, -i]) if a in bases]      #all bases at pos i
    if float(len(info))>=numsp*0.5:
        break
align = align[:, :-i]       #keep only bases up to the cutoff

for i in range(len(align[0])):      #same as above but from front
    info=[a for a in list(align[:, i]) if a in bases]
    if float(len(info))>=numsp*0.5:
        break
align = align[:, i:]

sp_to_remove=[]
for record in align:
    info = [p for p in record.seq if p in bases]
    if len(info) == 0:
        sp_to_remove.append(record.id)

#make new alignment wo missing seqs
if len(sp_to_remove) < len(align):
    newalign = MultipleSeqAlignment([record for record in align if record.id not in sp_to_remove], generic_dna)
    AlignIO.write(newalign, f+'2', fformat)