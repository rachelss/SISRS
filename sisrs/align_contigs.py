import os
import shutil
from glob import glob
from process import Process
from multiprocessing import Pool
from subprocess import Popen, PIPE, check_call


def sam_index_directory(dir_):
    basename = os.path.basename(dir_)
    filename = basename + '.bam'
    file_path = os.path.join(dir_, filename)

    index_command = [
        'samtools', 'index', file_path
    ]
    index_proc = Process(index_command)

    # TODO: should maybe be calling index_proc.wait() here...


class AlignContigsCommand(object):

    def __init__(self, data):

        self._data = data 

    def run(self):

        print("==== Renaming Scaffolds for SISRS ====")
        contig_dir = self._data['contig_dir']
        dir_lists = self._data['dir_lists']
        num_processors = self._data['num_processors']

        contig_file_path = os.path.join(contig_dir, 'contigs.fa')

        # backup contigs.fa
        backup_contig_file_path = os.path.join(
            contig_dir, 'contigs_OriginalNames.fa')
        shutil.move(contig_file_path, backup_contig_file_path)

        # rename scaffolds
        rename_command = [
            'rename.sh',
            'in={}'.format(backup_contig_file_path),
            'out={}'.format(contig_file_path),
            'prefix=SISRS',
            'addprefix=t'
        ]
        check_call(rename_command)
        print("==== Scaffolds Renamed ====")

        all_dirs = dir_lists.get_all_dirs()
        contig_prefix = os.path.join(contig_dir, 'contigs')

        # Build index
        build_command = [ 'bowtie2-build', contig_file_path, contig_prefix ]
        Process(build_command).wait()

        for dir_ in all_dirs:

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
            #os.remove(temp_file_path)
            #os.remove(header_file_path)

        print("==== Done Aligning ====")
        pool = Pool(num_processors)
        pool.map(sam_index_directory, all_dirs)
        print("==== Done Indexing Bam Files ====")


