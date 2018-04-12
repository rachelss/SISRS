# SISRS

SISRS: Site Identification from Short Read Sequences  
Version 1.6.2  
Copyright (c) 2013-2016 Rachel Schwartz <Rachel.Schwartz@asu.edu>  
https://github.com/rachelss/SISRS  
More information: Schwartz, R.S., K.M Harkins, A.C. Stone, and R.A. Cartwright. 2015. A composite genome approach to identify phylogenetically informative data from next-generation sequencing. BMC Bioinformatics. 16:193.
(http://www.biomedcentral.com/1471-2105/16/193/)

Talk from Evolution 2014 describing SISRS and its application:  
https://www.youtube.com/watch?v=0OMPuWc-J2E&list=UUq2cZF2DnfvIUVg4tyRH5Ng


# Requirements
* Built-In Genome Assemblers (Required if SISRS is building your composite genome)
  * Velvet (tested with v.1.2.10) (http://www.ebi.ac.uk/~zerbino/velvet/)
  * Minia (tested with v.2.0.7) (http://minia.genouest.org/)
  * AbySS (tested with v.2.0.2) (http://www.bcgsc.ca/platform/bioinfo/software/abyss)
* Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Python 3.6+, Biopython
* Samtools v1.3.1 (http://www.htslib.org/)
* MAFFT (http://mafft.cbrc.jp/alignment/software/)
* BBMap [requires Java] (https://sourceforge.net/projects/bbmap/)


# Running SISRS

## Install Docker
Follow the instructions [here](https://docs.docker.com/install/) to install
Docker **CE** for your operating system.

There's quite a bit going on in those instructions. As an example, if you're
on Ubuntu you'll follow the specific instructions
[here](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

Ignore anything that talks about Docker **EE**. They're just trying to sell you
stuff.

## Getting the Docker image

Download the SISRS docker image which comes with all the dependencies for
running SISRS.

`docker pull anderspitman/sisrs`

## Getting it up and running

First start up a docker container using the image obtained above:

`docker run -it anderspitman/sisrs bash`

Then from within the Docker container:

```
pip install sisrs
sisrs-python
```

## Running on example data

A (relatively) small example dataset is provided. This is the same data we
use for automated testing. First, retrieve the data by running the following
from within a running SISRS Docker container:

```
git clone --depth=1 https://github.com/anderspitman/SISRS_ci_data
```

Then try out some example commands (you'll need to remove the `output`
directory before running each command):

```
sisrs-python alignContigs -a premade -c 0 -f SISRS_ci_data/0_RawData_PremadeGenome -z output
sisrs-python identifyFixedSites -a premade -c 0 -f SISRS_ci_data/1_alignContigs -z output
sisrs-python outputAlignment -a premade -c 0 -f SISRS_ci_data/2_identifyFixedSites -z output
sisrs-python changeMissing -a premade -c 0 -f SISRS_ci_data/3_outputAlignment -z output
```


# Developing SISRS

The recommended way to develop SISRS is using the provided Docker image.
First, follow the instructions in the previous section to obtain the image.

Once you've downloaded the Docker image, you need to launch a container. The
recommended way to do this is to create the container with a mount point to
your SISRS directory. This allows you to make changes to the code from the host
system (your computer), while being able to run SISRS from within the Docker container. To
accomplish this, first navigate to the SISRS directory, then run the following
command:

```
docker run -it --mount src=$PWD,dst=/SISRS,type=bind -w /SISRS anderspitman/sisrs bash
```

This launches the container with your current directory mounted at /SISRS
within the container and starts bash inside that directory. Next, from within
the container, install SISRS in dev mode:

```
pip install -e .
```

This allows you to make changes to the code and run SISRS with the changes
immediately, without having to run pip install again.

You can verify your installation by running the test suite:

```
pytest -s tests
```


Input
=====

Next-gen sequence data such as Illumina HiSeq reads.
Data must be sorted into folders by taxon (e.g. species or genus).
Paired reads in fastq format must be specified by _R1 and _R2 in the (otherwise identical) filenames.
Paired and unpaired reads must have a fastq file extension.

Output
======

Nexus file with variable sites in a single alignment. Usable in most major phylogenetics software as a concatenated alignment with a setting for variable-sites-only.

Test Data
=========

The folder test_data (https://github.com/rachelss/SISRS_test_data) contains simulated data for 10 species on the tree found in simtree.tre . Using 40 processors this run took 9 minutes. Analysis of the alignment output by sisrs using raxml produced the correct tree.
