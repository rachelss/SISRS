SISRS
=====

SISRS: Site Identification from Short Read Sequences  
Version 1.6  
Copyright (c) 2013-2015 Rachel Schwartz <Rachel.Schwartz@asu.edu>  
https://github.com/rachelss/SISRS  
More information: http://arxiv.org/abs/1305.3665  

Talk from Evolution 2014 describing SISRS and its application:  
https://www.youtube.com/watch?v=0OMPuWc-J2E&list=UUq2cZF2DnfvIUVg4tyRH5Ng

License
=======

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

Requirements
============

* Velvet (including the perl script for merging paired reads) - http://www.ebi.ac.uk/~zerbino/velvet/
* Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* Python 2.7 and Biopython
* Samtools (http://www.htslib.org/)
* GNU Parallel (http://www.gnu.org/software/parallel/)

Input
=====

Next-gen sequence data such as Illumina HiSeq reads.
Data must be sorted into folders by taxon (e.g. species or genus).
Paired reads in fastq format must be specified by R1 and R2 in the (otherwise identical) filenames.
There should be no other R1 or R2 in the filename. If your file name is SRR287599_R2 you must change it to SRr287599_R2 or anything else you prefer.
Paired and unpaired reads must have a fastq file extension.

Running SISRS
=============

Usage:

 sisrs command options

By default, SISRS assumes that

 * A reference genome is not available.
 * The K-mer size to be used by Velvet in contig assembly is 21.
 * Only one processor is available.
 * Files are in fastq format.
 * A site is only required to have data for two species to be included
   in the final alignment.
 * Folders containing reads are in the present working directory.
 * A minimum of three reads are required to call the base at a site
   for a taxon.

Commands:

 * sites : produce an alignment of sites from raw reads
 * alignContigs : run sisrs skipping the composite genome assembly
 * mapContigs : run sisrs, also skipping alignment of reads to composite genome
 * identifyFixedSites : run sisrs, also skipping mapping of contigs to a reference
 * outputAlignment : get sisrs alignment from sites id'd for individual species
 * loci : produce a set of aligned loci based on the most variable regions of
        the composite genome
 
Flags:
    
 * -g : MANDATORY if running sisrs from the beginning - the approximate genome size
      - this will reduce the size of the composite assembly by using a subset
      of reads to approximate 10x coverage
 * -p : use this number of processors
 * -r : the path to the reference genome in fasta format
 * -k : k-mer size (for assembly)
 * -f : the folder containing the folders of reads
 * -n : the number of reads required to call a base at a site
 * -t : the threshold for calling a site; e.g. 0.99 means that >99% of
      bases for that taxon must be one allele; only recommended for
      low ploidy with <3 individuals
 * -m : the number of species that are allowed to have missing data at
      a site
 * -l : the number of alleles for sisrs loci
 * -a : change the assembler (specify Velvet or minia)

Output
======

Nexus file with variable sites in a single alignment. Usable in most major phylogenetics software as a concatenated alignment with a setting for variable-sites-only.

Test Data
=========

The folder test_data contains simulated data for 10 species on the tree found in simtree.tre . Using 40 processors this run took 9 minutes. Analysis of the alignment output by sisrs using raxml produced the correct tree.

Sample commands
==============

1. Basic sisrs run: start with fastq files and produce an alignment of variable sites 

   sisrs sites -g 1745690

2. Basic sisrs run with modifications

   sisrs sites -g 1745690 -p 40 -m 4 -f test_data -t .99 -a minia
 
2. Produce an alignment of loci based on the most variable loci in your basic sisrs run. Note - this command will run sisrs sites if (and only if) it was not run previously.

   sisrs loci -g 1745690 -p 40 -l 2 -f .
   
3. Get loci from your fastq files given known loci.

   first name your reference loci ref_genes.fa and put in your main folder
   
   sisrs loci -p 40 -f test_data
