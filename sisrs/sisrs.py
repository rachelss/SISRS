import subprocess
import os
import sys
import click
import shutil
from click import echo
from glob import glob
#import sub_sample_for_velvet_unshuff

def setup_output_directory(data_directory, output_directory, overwrite):

    if output_directory is None:
        output_directory = data_directory
        echo("Note: SISRS writing into data folder")
    elif not os.path.exists(output_directory):
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
@click.pass_context
def cli(ctx, assembler, data_directory, output_directory, overwrite):

    if data_directory is None:
        data_directory = os.getcwd()
    data_directory = os.path.abspath(data_directory)

    output_directory = setup_output_directory(
        data_directory, output_directory, overwrite)

    ctx.obj['data_dir'] = data_directory
    ctx.obj['out_dir'] = output_directory

@cli.command()
@click.option('--genome-size', '-g', required=True, type=int, help='Genome size')
@click.pass_context
def subsample(ctx, genome_size):
    #print(ctx.obj['data_dir'])
    pass

@cli.command()
@click.pass_context
def output_alignment(ctx):
    pass

def main():
    cli(obj={})
