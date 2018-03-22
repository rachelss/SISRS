import os
import shutil
from multiprocessing import Pool
from subprocess import check_call
from .process import Process 
from .aligners import create_aligner
from .command import Command


def sam_index_directory(dir_):
    basename = os.path.basename(dir_)
    filename = basename + '.bam'
    file_path = os.path.join(dir_, filename)

    index_command = [
        'samtools', 'index', file_path
    ]
    index_proc = Process(index_command)
    index_proc.wait()


class AlignContigsCommand(Command):

    def run(self):

        print("==== Renaming Scaffolds for SISRS ====")
        contig_dir = self._args['contig_dir']
        dir_lists = self._args['dir_lists']
        num_processors = self._args['num_processors']

        aligner = create_aligner(num_processors=num_processors)

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

        aligner.index(contig_file_path, contig_prefix)

        for dir_ in all_dirs:

            aligner.align(dir_, contig_prefix)

        print("==== Done Aligning ====")
        pool = Pool(num_processors)
        pool.map(sam_index_directory, all_dirs)
        print("==== Done Indexing Bam Files ====")


