#!/usr/local/bin/bash
#you must compile velvet prior to running sisrs - see the manual for appropriate compilation parameters
#run format: ./sisrs.sh plus flags
#folder should contain subfolders for each taxon of interest - subfolders should contain fastq (or fasta) files for that taxon
#reads must be paired, where the filename contains either R1 or R2

#flags
#-r <reference> #This is a file containing an assembled reference genome in fasta format
    #the first line of each sequence should only contain the sequence name
#-k <smallest kmer:largest kmer>    #optimize velvet for the best kmer size to produce the best reference contigs
    #default is kmer=21 - currently only uses default
#-f fasta or fastq
#-p <number> #number of processors that can be used
#-m <number> #number of species allowed to be missing data in the final alignment

#handle flags
while getopts r:v:p:f:m:a: option
do
    case "${option}"
    in
        r) REFFILE=${OPTARG};;
        k) KMER=${OPTARG};;
        p) PROCESSORS=${OPTARG};;
        f) FORMAT=${OPTARG};;
        m) MISSING=${OPTARG};;
        a) MAINFOLDER=${OPTARG};;
        esac
done
#required flags
assertNotEmpty() {
    : "${!1:? "$1 is empty, aborting."}"
}
assertNotEmpty FORMAT
assertNotEmpty MAINFOLDER
assertNotEmpty PROCESSORS
assertNotEmpty MISSING

MAINFOLDER=`echo "${MAINFOLDER}" | sed -e "s/\/*$//" `
VELVETLOC=$( dirname $( find / -name velvetg 2> /dev/null ) )
SHUFFLELOC=$( dirname $( find / -name shuffleSequences_fastq.pl 2> /dev/null ) )
FILELIST=$( find ${MAINFOLDER} -name *R1*${FORMAT} )
FOLDERLIST=$( dirname ${FILELIST} | sort -u | uniq )

#make reference contigs - velveth
#shuffle paired end reads together first
for FILE in $FILELIST; do
    ( perl ${SHUFFLELOC}/shuffleSequences_${FORMAT}.pl ${FILE} $(echo ${FILE}|sed 's/R1/R2/') $(echo ${FILE}|sed 's/R1/shuffled/') ) &
done
wait

if [ -z "$KMER" ]; then
    KMER=21
fi

if [ -n "${REFFILE}" ]; then
    ${VELVETLOC}/velveth ${MAINFOLDER}/velvetoutput ${KMER} -create_binary -${FORMAT} -shortPaired ${MAINFOLDER}/Sample*/*shuffled*${FORMAT} -fasta -long ${REFFILE}
else
    ${VELVETLOC}/velveth ${MAINFOLDER}/velvetoutput ${KMER} -create_binary -${FORMAT} -shortPaired ${MAINFOLDER}/Sample*/*shuffled*${FORMAT}
fi
echo ==== Velveth is finished ====

#run velvet - output contigs
${VELVETLOC}/velvetg ${MAINFOLDER}/velvetoutput 
echo ==== Velvetg is finished ====

#index contigs
bowtie2-build ${MAINFOLDER}/velvetoutput/contigs.fa ${MAINFOLDER}/velvetoutput/contigs

#align reads to contigs - run one processor per, in background
for FILE in $FILELIST; do
    NAME=$( echo ${FILE} | sed 's/R1//' | sed 's/\.[^.]*$//' )      #includes folder but not the read or the extension
    echo ==== Aligning ${NAME} ====
    #N=1 allows a mismatch  #x The basename of the index for the reference genome
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
    bwa index -a bwtsw ${REFFILE}    #index reference first
    bwa bwasw -t ${PROCESSORS} ${REFFILE} ${MAINFOLDER}/velvetoutput/contigs.fa > ${MAINFOLDER}/velvetoutput/align_contigs.sam
fi   
echo ==== Done Piling Up and Mapping Contigs ====

#put base for each site in a dictionary (allows no variation when calling sites)
for FOLDER in $FOLDERLIST; do
    ( python make_pruned_dict.py ${FOLDER} ) &    
done
wait
echo ==== Done Identifying Fixed Sites Without Error ====

if [ -n "${REFFILE}" ]; then
    ( python get_alignment.py ${MISSING} ${REFFILE} ${MAINFOLDER} ) &
else
    ( python get_alignment.py ${MISSING} X ${MAINFOLDER} ) &
fi
wait
echo ==== DONE ====
