import os
import sys
from multiprocessing import Pool

from .command import Command
from .sub_sample_for_velvet_unshuff import main as sub_sample_for_velvet_unshuff

def run_subsample_on_dir(args):
    dir_ = args[0]
    left_reads = args[1]
    sub_sample_for_velvet_unshuff(dir_, left_reads, '_R1', '_R2')

class SubsampleCommand(Command):

    def run(self):

        args = self._args

        genome_size = args['genome_size']
        if genome_size is None:
            print("Error: genome size required for subsampling")
            sys.exit(1)

        print("Subsampling data")
        out_dir = args['out_dir']
        dir_lists = args['dir_lists']
        all_dirs = dir_lists.get_all_dirs()
        num_processors = args['num_processors']

        subsample_dir = os.path.join(out_dir, 'subsamples')
        if os.path.exists(subsample_dir):
            os.remove(subsample_dir)
        os.makedirs(subsample_dir)

        numer = 10 * genome_size
        # assumes 100bp reads - add option in future
        denom = 100 * 2 * len(all_dirs)
        left_reads = numer / denom

        pool_args = [ (dir_, left_reads) for dir_ in all_dirs ]
        pool = Pool(num_processors)

        try:
            pool.map(run_subsample_on_dir, pool_args)
        except Exception as e:
            print("sub_sample_for_velvet_unshuff.py failed")
            print(e)
            sys.exit(1)

