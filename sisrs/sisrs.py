from multiprocessing import Pool
import os
import sys
import click
import argparse
from click import echo
from glob import glob
from pprint import pprint
from .align_contigs import AlignContigsCommand 
from .identify_fixed_sites import IdentifyFixedSitesCommand
from .get_alignment import main as get_alignment
from .change_missing import ChangeMissingCommand
from .process import Process


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
            sys.exit(1)

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

#@cli.command()
#@click.pass_context
def align_contigs(args_dict):

    command = AlignContigsCommand(args_dict)
    command.run()

@cli.command()
@click.option('--min-read', '-n', required=False, type=int, default=3,
              help='Min read')
@click.option('--threshold', '-t', required=False, type=int, default=1,
              help="Threshold for calling the site")
@click.pass_context
def identify_fixed_sites(ctx, min_read, threshold):

    ctx.obj['min_read'] = min_read
    ctx.obj['threshold'] = threshold 

    command = IdentifyFixedSitesCommand(ctx.obj)
    command.run()

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

@cli.command()
@click.pass_context
@click.option('--missing', '-m', required=False, type=int, 
              help="Num missing")
def change_missing(ctx, missing):

    command = ChangeMissingCommand(ctx.obj, missing)
    command.run()

def main():
    #cli(obj={})
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    align_contigs_parser = subparsers.add_parser(
            'alignContigs', help='Align Contigs')
    align_contigs_parser.add_argument('--data_directory')
    align_contigs_parser.add_argument('--output_directory')
    align_contigs_parser.add_argument('--overwrite')
    align_contigs_parser.add_argument(
            '-a', '--assembler', type=str, default='velvet')
    align_contigs_parser.add_argument('--num_processors', type=int, default=1)
    align_contigs_parser.set_defaults(func=align_contigs)

    change_missing_parser = subparsers.add_parser(
            'changeMissing', help='Change Missing')

    args = parser.parse_args()

    data_directory = args.data_directory
    output_directory = args.output_directory
    overwrite = args.overwrite
    assembler = args.assembler
    num_processors = args.num_processors

    if data_directory is None:
        data_directory = os.getcwd()
    data_directory = os.path.abspath(data_directory)

    output_directory = setup_output_directory(
        data_directory, output_directory, overwrite)

    dir_lists = DirectoryLists(output_directory)

    contig_dir = ''
    if assembler == 'premade':
        contig_dir = 'premadeoutput'

    args_dict = {}
    args_dict['data_dir'] = data_directory
    args_dict['out_dir'] = output_directory
    args_dict['assembler'] = assembler
    args_dict['dir_lists'] = dir_lists
    args_dict['contig_dir'] = os.path.join(output_directory, contig_dir)
    args_dict['num_processors'] = num_processors

    args.func(args_dict)
