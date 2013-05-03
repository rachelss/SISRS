
SISRS
=====

SISRS: SNP Identification from Short Read Sequences.

Copyright: Rachel Schwartz

License
=======

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

Requirements
============

* Velvet (including the perl script for merging paired reads)
* Bowtie2
* Python 2.7

Input
=====

Next-gen sequence data such as Illumina HiSeq reads.
Data (fasta or fastq files) should be sorted into folders by taxon (e.g. species or genus).
Data must be paired reads specified by R1 and R2 in the (otherwise identical) filenames.

Running SISRS
=============

SISRS can be run as 'bash sisrs.sh' from the command line. By default, SISRS assumes that
* A reference genome is not available.
* The K-mer size to be used by Velvet in contig assembly is 21.
* Only one processor is available.
* Files are in fastq format.
* A site is only required to have data for two species to be included in the final alignment.
* Folders containing reads are in the present working directory.
* A minimum of three reads are required to call the base at a site for a taxon.

Default settings can be changed using the following flags:
* -r : use to specify the location of the reference genome (must be a fasta file)
* -k : use to specify k-mer size
* -p : use to specify the number of processors available
* -f : use to specify reads in fasta format
* -m : use to specify the number of species allowed to have missing data at a site
* -a : use to specify the folder containing the folders of reads
* -n : use to specify the number of reads required to call a base at a site
* -s : use to skip all steps except the one identifying whether sites are variable among taxa

Example command: bash sisrs.sh -r ./reference.fasta -p 40 -f fastq -m 4 -a ./fastq_data/

Output
======

NEXUS file with SNPs in a single alignment. Usable in most major phylogenetics software.
