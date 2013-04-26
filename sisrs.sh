#!/usr/local/bin/bash
#you must compile velvet prior to running sisrs - see the manual for appropriate compilation parameters
#you must install bowtie2
#you must install samtools
#this script uses python2.7 - it has not been tested with python3
#run format: bash sisrs plus flags
#reads must be paired, where the filename contains either R1 or R2

#optional flags
#-r <reference> #This is a file containing an assembled reference genome in fasta format
    #the first line of each sequence should only contain the sequence name
#-k <kmer size> #optimize velvet for the best kmer size to produce the best reference contigs
    #default is kmer=21
#-p <number> #number of processors that can be used #default=1
#-f fasta or fastq #default=fastq
#-m <number> #number of species allowed to be missing data in the final alignment #default=total-2
#-a <folder> #folder should contain subfolders for each taxon of interest - subfolders should contain fastq/a files for that taxon #default = current dir
#-n <number> #number of reads required for a taxon to call that site - default is 3
#-s #will skip everything but comparing sites across taxa - use this to try different amounts of missing data

#handle flags
while getopts r:v:p:f:m:a:s option
do
case "${option}"
    in
        r) REFFILE=${OPTARG};;
        k) KMER=${OPTARG};;
        p) PROCESSORS=${OPTARG};;
        f) FORMAT=${OPTARG};;
        m) MISSING=${OPTARG};;
        a) MAINFOLDER=${OPTARG};;
        n) MINREAD=${OPTARG};;
        s) SKIP=1;;
        esac
done
#required flags
assertNotEmpty() {
    : "${!1:? "$1 is empty, aborting."}"
}

if [ -z "$KMER" ]; then
    KMER=21
fi

if [ -z "$MINREAD" ]; then
    MINREAD=3
fi

if [ -z "$PROCESSORS" ]; then
    PROCESSORS=1
fi

if [ -z "$FORMAT" ]; then
    FORMAT=fastq
fi

if [ -z "$MAINFOLDER" ]; then
    MAINFOLDER=.
fi

MAINFOLDER=`echo "${MAINFOLDER}" | sed -e "s/\/*$//" `
VELVETLOC=$( dirname $( find / -name velvetg 2> /dev/null ) )
SHUFFLELOC=$( dirname $( find / -name shuffleSequences_fastq.pl 2> /dev/null ) )
FILELIST=$( find ${MAINFOLDER} -name *R1*${FORMAT} )
FOLDERLIST=$( dirname ${FILELIST} | sort -u | uniq )

if [ -z "$MISSING" ]; then
    MISSING=0
    for i in $FOLDERLIST
    do
        MISSING=$(($MISSING + 1))
    done
    MISSING=$(($MISSING - 2))
fi

#use skip flag to skip everything but SNP finding - eg for different amounts of missing data
if [ -z "$SKIP" ]; then

#make reference contigs - velveth
#shuffle paired end reads together first
for FILE in $FILELIST; do
    ( perl ${SHUFFLELOC}/shuffleSequences_${FORMAT}.pl ${FILE} $(echo ${FILE}|sed 's/R1/R2/') $(echo ${FILE}|sed 's/R1/shuffled/') ) &
done
wait

if [ -n "${REFFILE}" ]; then
    ${VELVETLOC}/velveth ${MAINFOLDER}/velvetoutput ${KMER} -create_binary -fasta -reference ${REFFILE} -${FORMAT} -shortPaired ${MAINFOLDER}/*/*shuffled*${FORMAT}
else
    ${VELVETLOC}/velveth ${MAINFOLDER}/velvetoutput ${KMER} -create_binary -${FORMAT} -shortPaired ${MAINFOLDER}/*/*shuffled*${FORMAT}
fi
echo ==== Velveth is finished ====

#run velvet - output contigs
${VELVETLOC}/velvetg ${MAINFOLDER}/velvetoutput
echo ==== Velvetg is finished ====

#index contigs
bowtie2-build ${MAINFOLDER}/velvetoutput/contigs.fa ${MAINFOLDER}/velvetoutput/contigs

#align reads to contigs - run one processor per, in background
for FILE in $FILELIST; do
    NAME=$( echo ${FILE} | sed 's/R1//' | sed 's/\.[^.]*$//' ) #includes folder but not the read or the extension
    echo ==== Aligning ${NAME} ====
    #N=1 allows a mismatch #x The basename of the index for the reference genome
    ( bowtie2 -p 1 -N 1 --local -x ${MAINFOLDER}/velvetoutput/contigs -1 ${FILE} -2 $( echo ${FILE}|sed 's/R1/R2/' ) > >(tee ${NAME}_stdout.log) 2> >(tee ${NAME}_stderr.log >&2) | gzip > ${NAME}.sam.gz ) &
done
wait
echo ==== Done Aligning ====

#make bam files from zipped sam
for FILE in $FILELIST; do
NAME=$( echo ${FILE} | sed 's/R1//' | sed 's/\.[^.]*$//' )
    ( samtools view -uS ${NAME}.sam.gz | samtools sort -m 17179869184 - ${NAME};\
       samtools index ${NAME}.bam ) &
done
wait
echo ==== Done Making Bam Files ====

#get pileups for data
for FILE in $FILELIST; do
    NAME=$( echo ${FILE} | sed 's/R1//' | sed 's/\.[^.]*$//' )
    ( samtools mpileup -f ${MAINFOLDER}/velvetoutput/contigs.fa ${NAME}.bam > ${NAME}.pileups ) &
done
wait

#map contigs to reference
if [ -n "${REFFILE}" ]; then
    NAME=$( echo ${REFFILE} | sed 's/\.[^.]*$//' )
    bowtie2-build ${REFFILE} ${NAME}        #bowtie2-build [options]* <reference_in> <bt2_base>
    ( bowtie2 -p ${PROCESSORS} -N 1 -x ${NAME} -f -U ${MAINFOLDER}/velvetoutput/contigs.fa > >(tee ${NAME}_stdout.log) 2> >(tee ${NAME}_stderr.log >&2) -S ${MAINFOLDER}/velvetoutput/align_contigs.sam )    #bowtie2 -x <ref_base> -U <fq files> -S <output sam>
    bwa index -a bwtsw ${REFFILE} #index reference first
    bwa bwasw -t ${PROCESSORS} ${REFFILE} ${MAINFOLDER}/velvetoutput/contigs.fa > ${MAINFOLDER}/velvetoutput/align_contigs_bwa.sam
fi
echo ==== Done Piling Up and Mapping Contigs ====

#put base for each site in a dictionary (allows no variation when calling sites)
for FOLDER in $FOLDERLIST; do
    ( python get_pruned_dict_v1.py ${FOLDER} ${MINREAD} ) &
done
wait
echo ==== Done Identifying Fixed Sites Without Error ====

fi  #end of skip

if [ -n "${REFFILE}" ]; then
    ( python get_alignment_v1.py ${MISSING} ${REFFILE} ${MAINFOLDER} ) &
else
    ( python get_alignment_v1.py ${MISSING} X ${MAINFOLDER} ) &
fi
wait
echo ==== DONE ====
