import os
from glob import glob
from subprocess import PIPE 
from .aligner import Aligner
from ..process import Process

class Bowtie2Aligner(Aligner):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def index(self, reference_file_path, contig_path_prefix):
        build_command = [
            'bowtie2-build', reference_file_path, contig_path_prefix
        ]
        Process(build_command).wait()

    def align(self, directory, contig_prefix):

        dir_ = directory

        basename = os.path.basename(dir_)

        print("==== Aligning {} as Single-Ended ====".format(basename))
        fastq_filepaths = glob(os.path.join(dir_, '*.fastq'))

        # Generate temp file
        bowtie2_command = [
            'bowtie2',
            '-p', str(self._num_processors),
            '-N', '1',
            '--local',
            '-x', contig_prefix,
            '-U', ','.join(fastq_filepaths)
        ]
        bowtie2_proc = Process(bowtie2_command, stdout=PIPE)

        samtools_to_bam_command = [
            'samtools', 'view', '-Su',
            '-@', str(self._num_processors),
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
            '-@', str(self._num_processors),
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
            '-@', str(self._num_processors),
            '-H', temp_file_path,
            '-o', header_file_path
        ]
        samtools_header_proc = Process(samtools_header_command)

        samtools_header_proc.wait()


        # Generate final output file

        samtools_data_command = [
            'samtools', 'view',
            '-@', str(self._num_processors),
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
            '-@', str(self._num_processors),
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


