#!/usr/bin/env python2
from Bio import AlignIO
import sys
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

bases=['a','c','g','t','A','C','G','T']
f=sys.argv[1]
fformat = sys.argv[1].split('.')[-1]
align = AlignIO.read(f,fformat)

for i in range(1,len(align[0])):
    info=[a for a in list(align[:, -i]) if a in bases]
    if float(len(info))>float(len(align))*0.5:
        break
align = align[:, :-i]

for i in range(len(align[0])):
    info=[a for a in list(align[:, i]) if a in bases]
    if float(len(info))>float(len(align))*0.5:
        break
align = align[:, i:]

sp_to_remove=[]
for record in align:
    info = [p for p in record.seq if p in bases]
    if len(info) < 10:
        sp_to_remove.append(record.id)

#make new alignment wo short seqs
if len(sp_to_remove) < len(align):
    newalign = MultipleSeqAlignment([record for record in align if record.id not in sp_to_remove], generic_dna)
    AlignIO.write(newalign, f+'2', fformat)