import subprocess
from subprocess import Popen, PIPE
from multiprocessing import Pool
import os
import sys
import click
import shutil
from click import echo
from glob import glob
from get_alignment import main as get_alignment
from pprint import pprint
#import sub_sample_for_velvet_unshuff

class DirectoryLists(object):

    def __init__(self, base_dir):

        self._all_fastq = []
        
        for root, dirs, files in os.walk(base_dir):
            for filename in files:
                if filename.endswith('.fastq'):
                    file_path = os.path.join(root, filename)
                    self._all_fastq.append(file_path)

        self._all_dirs = sorted(list(set([ os.path.dirname(x) for x in self._all_fastq ])))

    def get_all_dirs(self):
        return self._all_dirs


def setup_output_directory(data_directory, output_directory, overwrite):

    if output_directory is None:
        output_directory = data_directory
        echo("Note: SISRS writing into data folder")
    elif not os.path.exists(output_directory):
        print(data_directory)
        recursive_symlinks(data_directory, output_directory)
    else:
        if overwrite:
            echo("{} already exists. Overwriting...".format(
                output_directory))
            recursive_symlinks(data_directory, output_directory)
        else:
            echo("{} already exists and overwrite flag not set. Aborting...".format(
                output_directory))

    return os.path.abspath(output_directory)

# based off https://stackoverflow.com/a/15824216/943814
def recursive_symlinks(src, dest):
    if os.path.isdir(src):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        files = glob(os.path.join(src, '*'))
        for f in files:
            f = os.path.basename(f)
            next_src = os.path.join(src, f)
            next_dest = os.path.join(dest, f)
            recursive_symlinks(
                next_src, next_dest)
    else:
        try:
            os.symlink(src, dest)
        except OSError as e:
            os.remove(dest)
            os.symlink(src, dest)


@click.group()
@click.option('--overwrite', type=bool, default=False)
@click.option(
    '--data-directory', '-f', type=click.Path(file_okay=False),
    help='Data directory')
@click.option(
    '--output-directory', '-z', type=click.Path(file_okay=False),
    help='Output directory')
@click.option('--assembler', '-a', default='velvet', help="Assembler")
@click.option('--continuous', '-c', default=1, help="Coninuous mode")
@click.option('--num-processors', '-p', default=1, help="Number of processors")
@click.pass_context
def cli(ctx, assembler, data_directory, output_directory, overwrite,
        continuous, num_processors):

    if data_directory is None:
        data_directory = os.getcwd()
    data_directory = os.path.abspath(data_directory)

    output_directory = setup_output_directory(
        data_directory, output_directory, overwrite)

    dir_lists = DirectoryLists(output_directory)

    contig_dir = ''
    if assembler == 'premade':
        contig_dir = 'premadeoutput'

    ctx.obj['data_dir'] = data_directory
    ctx.obj['out_dir'] = output_directory
    ctx.obj['assembler'] = assembler
    ctx.obj['dir_lists'] = dir_lists
    ctx.obj['contig_dir'] = os.path.join(output_directory, contig_dir)
    ctx.obj['num_processors'] = num_processors

@cli.command()
@click.option('--genome-size', '-g', required=True, type=int, help='Genome size')
@click.pass_context
def subsample(ctx, genome_size):
    #print(ctx.obj['data_dir'])
    pass

@cli.command()
@click.pass_context
def output_alignment(ctx):
    data_dir = ctx.obj['data_dir']
    out_dir = ctx.obj['out_dir']
    assembler = ctx.obj['assembler']
    dir_lists = ctx.obj['dir_lists']

    all_dirs = dir_lists.get_all_dirs()
    num_missing = len(all_dirs) - 2
    get_alignment(num_missing, 'X', out_dir, assembler)

def sam_index_directory(dir_):
    basename = os.path.basename(dir_)
    filename = basename + '.bam'
    file_path = os.path.join(dir_, filename)

    index_command = [
        'samtools', 'index', file_path
    ]
    index_proc = Popen(index_command)

@cli.command()
@click.pass_context
def align_contigs(ctx):

    print("==== Renaming Scaffolds for SISRS ====")
    contig_dir = ctx.obj['contig_dir']
    dir_lists = ctx.obj['dir_lists']
    num_processors = ctx.obj['num_processors']

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
    subprocess.check_call(rename_command)
    print("==== Scaffolds Renamed ====")

    all_dirs = dir_lists.get_all_dirs()
    contig_prefix = os.path.join(contig_dir, 'contigs')

    # Build index
    build_command = [ 'bowtie2-build', contig_file_path, contig_prefix ]
    subprocess.check_call(build_command)

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
        bowtie2_proc = Popen(bowtie2_command, stdout=PIPE)

        samtools_to_bam_command = [
            'samtools', 'view', '-Su',
            '-@', str(num_processors),
            '-F' , '4',
            '-',
        ]
        samtools_to_bam_proc = Popen(
            samtools_to_bam_command, stdin=bowtie2_proc.stdout,
            stdout=PIPE)

        output_base_path = os.path.join(dir_, basename)
        temp_file_path = output_base_path + '_Temp.bam'

        samtools_sort_command = [
            'samtools', 'sort',
            '-@', str(num_processors),
            '-',
            '-o', temp_file_path,
        ]
        samtools_sort_proc = Popen(
            samtools_sort_command, stdin=samtools_to_bam_proc.stdout)

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
        samtools_header_proc = Popen(samtools_header_command)

        samtools_header_proc.wait()


        # Generate final output file

        samtools_data_command = [
            'samtools', 'view',
            '-@', str(num_processors),
            temp_file_path,
        ]
        samtools_data_proc = Popen(samtools_data_command, stdout=PIPE)

        grep_command = [
            'grep', '-v', 'XS:'
        ]
        grep_proc = Popen(
            grep_command, stdin=samtools_data_proc.stdout, stdout=PIPE)

        cat_command = [
            'cat', header_file_path, '-'
        ]
        cat_proc = Popen(cat_command, stdin=grep_proc.stdout, stdout=PIPE)

        samtools_final_command = [
            'samtools', 'view',
            '-@', str(num_processors),
            '-b',
            '-',
            '-o', output_base_path + '.bam'
        ]

        samtools_final_proc = Popen(
            samtools_final_command, stdin=cat_proc.stdout)

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


def main():
    cli(obj={})