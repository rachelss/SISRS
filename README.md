SISRS
=====

SISRS: SNP Identification from Short Read Sequences.

Copyright: Rachel Schwartz

More information: http://arxiv.org/abs/1305.3665

Talk from Evolution 2014 describing SISRS and its application: https://www.youtube.com/watch?v=0OMPuWc-J2E&list=UUq2cZF2DnfvIUVg4tyRH5Ng

License
=======

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

Requirements
============

* Velvet (including the perl script for merging paired reads) - http://www.ebi.ac.uk/~zerbino/velvet/
* Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* Python 2.7 and Biopython
* Samtools

Input
=====

Next-gen sequence data such as Illumina HiSeq reads.
Data (fastq files) should be sorted into folders by taxon (e.g. species or genus).
Data must be paired reads specified by R1 and R2 in the (otherwise identical) filenames.

Running SISRS
=============

SISRS can be run as 'sisrs' from the command line if it is installed in a location listed in your path. Otherwise run SISRS as <directory>/sisrs. 

By default, SISRS assumes that
* A reference genome is not available.
* The K-mer size to be used by Velvet in contig assembly is 21.
* Only one processor is available.
* Files are in fastq format.
* A site is only required to have data for two species to be included in the final alignment.
* Folders containing reads are in the present working directory.
* A minimum of three reads are required to call the base at a site for a taxon.

Default settings can be changed using the following flags:
* -g : use to specify the approximate genome size - this is optional but will reduce the size of the composite assembly by using a subset of reads to approximate 10x coverage
* -r : use to specify the location of the reference genome (must be a fasta file)
* -k : use to specify k-mer size
* -p : use to specify the number of processors available
* -f : use to specify reads in fasta format
* -m : use to specify the number of species allowed to have missing data at a site
* -a : use to specify the folder containing the folders of reads
* -n : use to specify the number of reads required to call a base at a site
* -t : use to specify the threshold for calling a site; e.g. 0.99 means that >99% of bases for that taxon must be one allele; only recommended for low ploidy with <3 individuals
* -s : use to skip all steps except the one identifying whether sites are variable among taxa

Example command: bash sisrs.sh -r ./reference.fasta -p 40 -f fastq -m 4 -a ./fastq_data/

Output
======

NEXUS file with SNPs in a single alignment. Usable in most major phylogenetics software.
