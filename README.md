SISRS_B
=====

SISRS: Site Identification from Short Read Sequences  
Version 1.6  
Copyright (c) 2013-2016 Rachel Schwartz <Rachel.Schwartz@asu.edu>  
https://github.com/rachelss/SISRS  
More information: Schwartz, R.S., K.M Harkins, A.C. Stone, and R.A. Cartwright. 2015. A composite genome approach to identify phylogenetically informative data from next-generation sequencing. BMC Bioinformatics. 16:193.
(http://www.biomedcentral.com/1471-2105/16/193/)

Talk from Evolution 2014 describing SISRS and its application:  
https://www.youtube.com/watch?v=0OMPuWc-J2E&list=UUq2cZF2DnfvIUVg4tyRH5Ng

License
=======

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

Requirements
============
* Built-In Genome Assemblers (Required if SISRS is building your composite genome)
  * Velvet (tested with v.1.2.10) (http://www.ebi.ac.uk/~zerbino/velvet/)
  * Minia (tested with v.2.0.7) (http://minia.genouest.org/)
  * AbySS (tested with v.2.0.2) (http://www.bcgsc.ca/platform/bioinfo/software/abyss)
* Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Python 2.7, Biopython, and PySAM
* Samtools v1.3.1 (http://www.htslib.org/)
* GNU Parallel (http://www.gnu.org/software/parallel/)
* MAFFT (http://mafft.cbrc.jp/alignment/software/)
* BBMap [requires Java] (https://sourceforge.net/projects/bbmap/)


Input
=====

Next-gen sequence data such as Illumina HiSeq reads.
Data must be sorted into folders by taxon (e.g. species or genus).
Paired reads in fastq format must be specified by _R1 and _R2 in the (otherwise identical) filenames.
Paired and unpaired reads must have a fastq file extension.

Running SISRS
=============

## Usage:

 sisrs command options

 ### By default, SISRS assumes that

  * A reference genome is not available and a composite assembly will be
    assembled using Velvet
  * The K-mer size to be used by Velvet in contig assembly is 21.
  * Only one processor is available.
  * Files are in fastq format.
  * Paired read filenames end with _R1 and _R2
  * A site is only required to have data for two species to be included in the
    final alignment.
  * Folders containing reads are in the present working directory
  * SISRS data will be output into the present working directory
  * A minimum of three reads are required to call the base at a site
    for a taxon.

### Commands:  

**sites**: produce an alignment of sites from raw reads  

**loci**: produce a set of aligned loci based on the most variable regions of the composite genome  


#### Subcommands of sites:

**subSample**: run sisrs subsampling scheme, subsampling reads from all taxa to ~10X coverage across species, relative to user-specified genome size  

**buildContigs**: given subsampled reads, run sisrs composite genome assembly with user-specified assembler  

**alignContigs**: align reads to composite genome as single-ended, uniqely mapped  

**mapContigs**: align composite genome reads to a reference genome (optional)  

**identifyFixedSites**: find sites with no within-taxa variation  

**outputAlignment**: output alignment file of sisrs sites  

**changeMissing**: given alignment of sites (alignment.nex), output a file with only sites missing fewer than a specified number of samples per site  


 #### Option Flags:

 * -g : the approximate genome size (**mandatory** if sisrs will be assembling a composite genome)
 * -p : use this number of processors *[Default: 1]*
 * -r : the path to the reference genome in fasta format *[Optional]*
 * -k : k-mer size (for assembly) *[Default: 21]*  
 * -f : absolute path to the directory containing the folders of reads *[Default: Current Directory]*
 * -z : absolute path to either empty or non-existent directory where SISRS will output data *[Default: Current Directory]*
 * -n : the number of reads required to call a base at a site  *[Default: 3]*
 * -t : the threshold for calling a site; e.g. 0.99 means that >99% of bases for that taxon must be one allele; only recommended for low ploidy with <3 individuals  *[Default: 1 (100%)]*
 * -m : the number of species that are allowed to have missing data at a site
 * -o : the length of the final loci dataset for dating  
 * -l : the number of alleles  
 * -a : assembler [velvet, minia, abyss, or premade; *Default: velvet*]
      - If using a premade composite genome, it must be in a folder named 'premadeoutput' in the same directory as the folders of read data, and must be called 'contigs.fa'
 * -s : Sites to analyze when running 'loci' [0,1,2]
      - 0 [Default], all variable sites, including singletons
      - 1, variable sites excluding singletons
      - 2, only biallelic variable sites
 * -c : continuous command mode for calling subcommands [1,0]  
      - 1 [Default]: calling a subcommand runs that subcommand **and all subsequent steps in the pipeline**
      - 0: calling a subcommand runs **only** that subcommand

Output
======

Nexus file with variable sites in a single alignment. Usable in most major phylogenetics software as a concatenated alignment with a setting for variable-sites-only.

Test Data
=========

The folder test_data (https://github.com/rachelss/SISRS_test_data) contains simulated data for 10 species on the tree found in simtree.tre . Using 40 processors this run took 9 minutes. Analysis of the alignment output by sisrs using raxml produced the correct tree.

Sample commands
==============

1. Basic sisrs run: start with fastq files and produce an alignment of variable sites
```
sisrs sites -g 1745690
```
2. Basic sisrs run with modifications
```
sisrs sites -g 1745690 -p 40 -m 4 -f /usr/test_data -z /usr/output_data -t .99 -a minia
```
3. Run only sisrs read subsampling step
```
sisrs subSample -g 1745690 -f /usr/test_data -c 0
```
4. Produce an alignment of loci based on the most variable loci in your basic sisrs run. Note - this command will run sisrs sites if (and only if) it was not run previously.
```
sisrs loci -g 1745690 -p 40 -l 2 -f /usr/test_data           # Will run sites first, then loci
sisrs loci -g 1745690 -p 40 -l 2 -f /usr/SISRS_sites_ouput   # Will run loci from previous sites data
```
5. Get loci from your fastq files given known loci.

   first name your reference loci ref_genes.fa and put in your main folder
```
sisrs loci -p 40 -f /usr/test_data
```
