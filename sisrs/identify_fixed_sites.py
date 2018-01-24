import os
import sys
from process import Process, AlignmentProcess
from multiprocessing import Pool
from specific_genome import main as specific_genome
from get_pruned_dict import main as get_pruned_dict

def run_mpileup(args):

    dir_ = args[0]
    contig_file_path = args[1]
    
    taxon_name = os.path.basename(dir_)

    path = os.path.join(dir_, taxon_name)

    print(path)

    command = [
        'samtools', 'mpileup',
        '-f', contig_file_path,
        '-o', path + '.pileups',
        path + '.bam', 
    ]
    Process(command).wait()

def run_specific_genome(args):

    dir_ = args[0]
    contig_file_path = args[1]

    specific_genome(dir_, contig_file_path)

def run_faidx(dir_):

    path = os.path.join(dir_, 'contigs.fa')

    command = [
        'samtools', 'faidx', path
    ]
    Process(command).wait()

def run_bowtie2(dir_):

    contigs_path = os.path.join(dir_, 'contigs.fa')
    contigs_dir = os.path.join(dir_, 'contigs')

    command = [
        'bowtie2-build', contigs_path, contigs_dir,
    ]
    Process(command).wait()

def run_index(dir_):

    taxon_name = os.path.basename(dir_)
    bam_path = os.path.join(dir_, taxon_name + '.bam')

    command = [
        'samtools', 'index',
        bam_path,
    ]
    Process(command).wait()

def run_get_pruned_dict(args):
    get_pruned_dict(*args)


class IdentifyFixedSitesCommand(object):

    def __init__(self, data):

        self._data = data 

    def run(self):

        data = self._data

        contig_dir = data['contig_dir']
        num_processors = data['num_processors']
        dir_lists = data['dir_lists']
        min_read = data['min_read']
        threshold = data['threshold']

        all_dirs = dir_lists.get_all_dirs()

        contig_file_path = os.path.join(contig_dir, 'contigs.fa')

        samtools_faidx_command = [
            'samtools', 'faidx', contig_file_path
        ]
        Process(['samtools', 'faidx', contig_file_path]).wait()


        pool = Pool(num_processors)
        args = [ (dir_, contig_file_path) for dir_ in all_dirs ]

        pool.map(run_mpileup, args)

        try:
            pool.map(run_specific_genome, args)
        except Exception:
            print("specific_genome.py failed")
            sys.exit(1)

        pool.map(run_faidx, all_dirs)
        pool.map(run_bowtie2, all_dirs)


        for dir_ in all_dirs:

            taxon_name = os.path.basename(dir_)
            os.remove(os.path.join(dir_, taxon_name + '.bam'))
            os.remove(os.path.join(dir_, taxon_name + '.bam.bai'))
            os.remove(os.path.join(dir_, taxon_name + '.pileups'))

            contig_prefix = os.path.join(dir_, 'contigs')
            AlignmentProcess(dir_, contig_prefix, num_processors)

        args = []
        for dir_ in all_dirs:
            contig_file_path = os.path.join(dir_, 'contigs.fa')
            args.append([dir_, contig_file_path])

        pool.map(run_index, all_dirs)
        pool.map(run_mpileup, args)
 
        # put base for each site in a dictionary (allows no variation when
        # calling sites)
        args = [ (dir_, min_read, threshold) for dir_ in all_dirs ]
        try:
            pool.map(run_get_pruned_dict, args)
        except Exception:
            print("get_pruned_dict.py failed")
            sys.exit(1)
        print("==== Done Identifying Fixed Sites Without Error ====")