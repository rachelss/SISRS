<<<<<<< HEAD
#!/usr/bin/env python
=======
#!/usr/bin/env python2
>>>>>>> release/1.6
import sys
from Bio import SeqIO

# Get arguments, required and optional
argc = len(sys.argv)
if argc < 3 or argc > 4:
    sys.exit('Usage: python get_seq.py CONTIG.fa LOCI.txt <output.fa>')
contig_file = sys.argv[1]
loci_file = sys.argv[2]
out_file = "ref_genes.fa"
if argc == 4:
    out_file = sys.argv[3]

# Read contig fasta file into dictionary with sequence ID as the key
fasta_dict = {}
contig_handle = open(contig_file, "r")
fasta_seq = SeqIO.parse(contig_handle, 'fasta')
for read in fasta_seq:
    name, seq = read.id, read.seq
    fasta_dict[name] = read
contig_handle.close()

# Go through each loci (taken from alignment.nex) and collect the corresponding
# read from the contig dictionary
ref_seq = []
loci_handle = open(loci_file, "r")
for line in loci_handle:
    ref_seq.append(fasta_dict[line.strip('\n')])

# write all matches between loci.txt and config.fa into new fasta file
out_handle = open(out_file, "w")
SeqIO.write(ref_seq, out_handle, "fasta")
out_handle.close()
