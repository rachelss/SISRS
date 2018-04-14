import os
import shutil
from glob import glob

from .assembler import Assembler
from ..process import Process


class MiniaAssembler(Assembler):

    def assemble(self):
        """Assember minia"""

        # NOTE: This is untested and incomplete. I think my version of minia
        # is a lot newer than the one being used in the bash script I ported
        # this from. At the very least the output files aren't being put in
        # the right places.

        self._subsample_dir = os.path.join(self._out_dir, 'subsamples')

        minia_out_dir = os.path.join(self._out_dir, 'miniaoutput')

        if not os.path.exists(minia_out_dir):
            os.mkdir(minia_out_dir)

        fastq_paths = glob(os.path.join(self._subsample_dir, "*fastq"))

        minia_read_file_path = os.path.join(self._subsample_dir, 'minia_read_file.txt')

        with open(minia_read_file_path, 'w') as f:
            f.write('\n'.join(fastq_paths))

        minia_command = [
            'minia',
            '-in', minia_read_file_path,
            '-kmer-size', str(self._kmer_size),
            '-nb-cores', str(self._num_processors),
            '-out-dir', minia_out_dir,
        ]

        Process(minia_command).wait()

        contigs_filename = 'minia_read_file.contigs.fa'
        contigs_path = os.path.join(os.getcwd(), contigs_filename)
        dest_contigs_path = os.path.join(minia_out_dir, 'contigs.fa')
        #shutil.move(contigs_path, dest_contigs_path)
        shutil.move(contigs_path, dest_contigs_path)

        temp_files = [
            os.path.join(minia_out_dir, 'minia_read_file.h5'),
        ] + glob("*unitigs*")

        for path in temp_files:
            os.remove(path)
