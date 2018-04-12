#!/usr/bin/env python3
import os
import sys
import shutil
from subprocess import check_call
from Bio import SeqIO

contig_dir = sys.argv[1]
contig_file_path = os.path.join(contig_dir, 'contigs.fa')

# backup contigs.fa
backup_contig_file_path = os.path.join(
    contig_dir, 'contigs_OriginalNames.fa')
shutil.move(contig_file_path, backup_contig_file_path)

# rename scaffolds
rename_command = [
    'rename.sh',
    'in={}'.format(backup_contig_file_path),
    'out={}'.format(contig_file_path),
    'prefix=SISRS',
    'addprefix=t'
]
check_call(rename_command)
print("==== Scaffolds Renamed ====",flush=True)

#CREATE FILE OF ALL CONTIG SEQUENCE LENGTHS

seqLengthFile = open(contig_dir+'/contigs_SeqLength.tsv', "w")
for seq_record in SeqIO.parse(contig_file_path,"fasta"):
  seqLengthFile.write(str(seq_record.id)+"\t"+str(len(seq_record))+"\n")
seqLengthFile.close()
print("==== Congig Length File Generated ====",flush=True)

#CREATE FILE WITH EVERY SITE IN ALIGNMENT
siteCount=0
locListFile = open(contig_dir+'/contigs_LocList','a+')
with open(contig_dir +"/contigs_SeqLength.tsv","r") as filein:
    for line in iter(filein):
        splitline=line.split()
        for x in range(1,(int(splitline[1])+1)):
            locListFile.write((splitline[0] +'/'+str(x)+'\n'))
            siteCount+=1
locListFile.close()
print("==== Site list created: " + str(siteCount) + " total sites ==== \n",flush=True)
