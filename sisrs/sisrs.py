from multiprocessing import Pool
import os
import sys
import argparse
from glob import glob
from pprint import pprint
from .subsample import SubsampleCommand
from .align_contigs import AlignContigsCommand 
from .identify_fixed_sites import IdentifyFixedSitesCommand
from .output_alignment import OutputAlignmentCommand
from .change_missing import ChangeMissingCommand
from .process import Process


class SISRSPipeline(object):
    def __init__(self):
         
        self._pipeline = [
            [ 'subSample', self._subsample ],
            [ 'alignContigs', self._align_contigs ],
            [ 'identifyFixedSites', self._identify_fixed_sites ],
            [ 'outputAlignment', self._output_alignment],
            [ 'changeMissing', self._change_missing],
        ]

        self._pipeline_lookup = {
            stage[0]: {
                'index': i,
                'func': stage[1]
            } for i, stage in enumerate(self._pipeline)
        }

    def get_command_list(self):
        commands = []
        for stage in self._pipeline:
            commands.append(stage[0])
        return commands

    def run_command(self, command_name, args):
        self._pipeline_lookup[command_name]['func'](args)

    def run_commands_from(self, first_command_name, args):
        
        first_command_index = \
            self._pipeline_lookup[first_command_name]['index']

        for stage in self._pipeline[first_command_index:]:
            command_name = stage[0]
            self.run_command(command_name, args)

    @staticmethod
    def _subsample(args):
        command = SubsampleCommand(args)
        command.run()

    @staticmethod
    def _align_contigs(args):
        command = AlignContigsCommand(args)
        command.run()

    @staticmethod
    def _identify_fixed_sites(args):
        command = IdentifyFixedSitesCommand(args)
        command.run()
        
    @staticmethod
    def _output_alignment(args):
        command = OutputAlignmentCommand(args)
        command.run()

    @staticmethod
    def _change_missing(args):
        command = ChangeMissingCommand(args)
        command.run()


class DirectoryLists(object):

    def __init__(self, base_dir):

        self._all_fastq = []
        self._paired = []
        self._unpaired = []
        
        for root, dirs, files in os.walk(base_dir):
            for filename in files:

                file_path = os.path.join(root, filename)

                if filename.endswith('.fastq'):
                    self._all_fastq.append(file_path)

                if self._is_paired_read_filename(filename):
                    self._paired.append(file_path)
                else:
                    self._unpaired.append(file_path)

        self._all_dirs = sorted(list(set([ os.path.dirname(x) for x in self._all_fastq ])))

    def get_all_dirs(self):
        return self._all_dirs

    def get_paired(self):
        return self._paired

    def get_unpaired(self):
        return self._unpaired

    def _is_paired_read_filename(self, filename):
        return (filename.endswith('_R1.fastq') or
                filename.endswith('_R2.fastq') or
                # TODO: is it true that subsampled files are always paired?
                'subsampled' in filename)


def setup_output_directory(data_directory, output_directory, overwrite):

    if output_directory is None:
        output_directory = data_directory
        print("Note: SISRS writing into data folder")
    elif not os.path.exists(output_directory):
        print(data_directory)
        recursive_symlinks(data_directory, output_directory)
    else:
        if overwrite:
            print("{} already exists. Overwriting...".format(
                output_directory))
            recursive_symlinks(data_directory, output_directory)
        else:
            print("{} already exists and overwrite flag not set. Aborting...".format(
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


def main():

    pipeline = SISRSPipeline()

    parser = argparse.ArgumentParser()

    parser.add_argument(
            'command',
            choices=pipeline.get_command_list(),
            help="Command to run")
    parser.add_argument('-c', '--continuous', default=1)
    parser.add_argument('-f', '--data-directory')
    parser.add_argument('-z', '--output-directory')
    parser.add_argument('--overwrite')
    parser.add_argument('-a', '--assembler', type=str, default='velvet')
    parser.add_argument('-p', '--num_processors', type=int, default=1)
    parser.add_argument('-m', '--missing', type=int, help="Num missing")
    parser.add_argument('-g', '--genome-size', type=int, help="Genome size",
            default=None)

    # identifyFixedSites specific
    parser.add_argument(
            '-n', '--min-read', required=False, type=int, default=3,
            help='Min read')
    parser.add_argument(
            '-t', '--threshold', required=False, type=int, default=1,
              help="Threshold for calling the site")

    args = parser.parse_args()

    data_directory = args.data_directory
    output_directory = args.output_directory
    overwrite = args.overwrite
    assembler = args.assembler
    num_processors = args.num_processors
    continuous = args.continuous

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
    args_dict['continuous'] = continuous
    args_dict['min_read'] = args.min_read 
    args_dict['threshold'] = args.threshold
    args_dict['missing'] = args.missing
    args_dict['genome_size'] = args.genome_size

    command_name = args.command

    if continuous == 1:
        pipeline.run_commands_from(command_name, args_dict)
    else:
        pipeline.run_command(command_name, args_dict)

