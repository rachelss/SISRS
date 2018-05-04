import os
from glob import glob

from .assembler import Assembler
from ..process import Process


class VelvetAssembler(Assembler):

    # NOTE: This is currently untested

    def __init__(self, dir_lists, out_dir, **kwargs):
        super().__init__(dir_lists, out_dir, **kwargs)

        self._ref_file_path = kwargs['ref_file_path']

    def assemble(self):
        """Assember velvet"""

        paired_paths = self._dir_lists.get_paired()
        len_paired = len(paired_paths)
        unpaired_paths = self._dir_lists.get_unpaired()
        len_unpaired = len(unpaired_paths)
        
        command_builder = VelvetCommandBuilder(self._out_dir, self._kmer_size)

        if self._ref_file_path is not None:

            command_builder.add_reference(self._ref_file_path)

            if len_unpaired > 0:
                if len_paired > 0:
                    print("Running Velvet with PE and SE reads, and reference")
                    command_builder.add_paired_end()
                    command_builder.add_single_end()
                else:
                    print("Running Velvet with SE reads and reference")
                    command_builder.add_single_end()
            else:
                print("Running Velvet with PE reads and reference")
        else:
            if len_unpaired > 0:
                if len_paired > 0:
                    print("Running Velvet with PE and SE reads")
                    command_builder.add_paired_end()
                    command_builder.add_single_end()
                else:
                    print("Running Velvet with SE reads")
                    command_builder.add_single_end()
            else:
                print("Running Velvet with PE reads")
                command_builder.add_paired_end()

        command = command_builder.build()
        Process(command).wait() 

        velvet_out_dir = os.path.join(self._out_dir, 'velvetoutput')
        velvetg_command = [
            'velvetg', velvet_out_dir,
            '-exp_cov', 'auto' '-cov_cutoff' 'auto'
        ]

        Process(velvetg_command).wait()


class VelvetCommandBuilder(object):

    def __init__(self, out_dir, kmer_size):

        velvet_out_dir = os.path.join(out_dir, 'velvetoutput')
        self._subsample_dir = os.path.join(out_dir, 'subsamples')

        self._command = [
            'velveth',
            velvet_out_dir,
            str(kmer_size),
            '-create_binary',
        ]

    def add_reference(self, ref_file_path):
        self._command += [ '-fasta', '-reference', ref_file_path ]
        return self

    def add_paired_end(self):

        paths = glob(os.path.join(self._subsample_dir, "*subsampledp.fastq"))
        paired_option = [ '-fastq', '-shortPaired' ] + paths
        self._command += paired_option
        return self

    def add_single_end(self):
        paths = glob(os.path.join(self._subsample_dir, "*subsampledu.fastq"))
        paired_option = [ '-fastq', '-short' ] + paths
        self._command += paired_option
        return self

    def build(self):
        return self._command
