import os
from subprocess import Popen
from glob import glob
from subprocess import Popen, PIPE 


class Process(object):

    def __init__(self, command, stdin=None, stdout=None):
        try:
            self.proc = Popen(command, stdin=stdin, stdout=stdout)
        except Exception as e:
            print("Error calling: {}".format(command[0]))
            raise e 

    def wait(self):
        return self.proc.wait()

    def pipe(self):
        return self.proc.stdout


class AlignmentProcess(object):

    def __init__(self, dir_, contig_prefix, num_processors):

        basename = os.path.basename(dir_)

        print("==== Aligning {} as Single-Ended ====".format(basename))
        fastq_filepaths = glob(os.path.join(dir_, '*.fastq'))

        # Generate temp file
        bowtie2_command = [
            'bowtie2',
            '-p', str(num_processors),
            '-N', '1',
            '--local',
            '-x', contig_prefix,
            '-U', ','.join(fastq_filepaths)
        ]
        bowtie2_proc = Process(bowtie2_command, stdout=PIPE)

        samtools_to_bam_command = [
            'samtools', 'view', '-Su',
            '-@', str(num_processors),
            '-F' , '4',
            '-',
        ]
        samtools_to_bam_proc = Process(
            samtools_to_bam_command, stdin=bowtie2_proc.pipe(),
            stdout=PIPE)

        output_base_path = os.path.join(dir_, basename)
        temp_file_path = output_base_path + '_Temp.bam'

        samtools_sort_command = [
            'samtools', 'sort',
            '-@', str(num_processors),
            '-',
            '-o', temp_file_path,
        ]
        samtools_sort_proc = Process(
            samtools_sort_command, stdin=samtools_to_bam_proc.pipe())

        bowtie2_proc.wait()
        samtools_to_bam_proc.wait()
        samtools_sort_proc.wait()


        # Generate header file

        header_file_path = output_base_path + '_Header.sam'

        samtools_header_command = [
            'samtools', 'view',
            '-@', str(num_processors),
            '-H', temp_file_path,
            '-o', header_file_path
        ]
        samtools_header_proc = Process(samtools_header_command)

        samtools_header_proc.wait()


        # Generate final output file

        samtools_data_command = [
            'samtools', 'view',
            '-@', str(num_processors),
            temp_file_path,
        ]
        samtools_data_proc = Process(samtools_data_command, stdout=PIPE)

        grep_command = [
            'grep', '-v', 'XS:'
        ]
        grep_proc = Process(
            grep_command, stdin=samtools_data_proc.pipe(), stdout=PIPE)

        cat_command = [
            'cat', header_file_path, '-'
        ]
        cat_proc = Process(cat_command, stdin=grep_proc.pipe(), stdout=PIPE)

        samtools_final_command = [
            'samtools', 'view',
            '-@', str(num_processors),
            '-b',
            '-',
            '-o', output_base_path + '.bam'
        ]

        samtools_final_proc = Process(
            samtools_final_command, stdin=cat_proc.pipe())

        samtools_data_proc.wait()
        grep_proc.wait()
        cat_proc.wait()
        samtools_final_proc.wait()

        # Remove temp files
        os.remove(temp_file_path)
        os.remove(header_file_path)


class MpileupProcess(object):

    def __init__(self, dir_, contig_file_path, num_processors):

        taxon_name = os.path.basename(dir_)
        bam_path = os.path.join(dir_, '{}.bam'.format(taxon_name))
        pileups_path  = os.path.join(dir_, '{}.pileups'.format(taxon_name))

        command = [
            'samtools',
            '-f', contig_file_path,
            bam_path,
            '-o', pileups_path
        ]
        Process(command).wait()

    #parallel --jobs "${PROCESSORS}" 'samtools mpileup -f' "${OUTFOLDER}"/"${CONTIGS}"/contigs.fa '"$( echo {}/$(basename {} ) )".bam' '> "$( echo {}/$(basename {} ) )".pileups' ::: "${FOLDERLISTA[@]}"
