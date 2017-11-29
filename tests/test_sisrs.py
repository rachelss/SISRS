import os
import subprocess
import shutil
from filecmp import cmp
from os.path import join, exists

data_base_dir = 'pipeline_stages'
out_base_dir = 'output'

print("out:")
print(out_base_dir)
if os.path.exists(out_base_dir):
    print("remove it")
    shutil.rmtree(out_base_dir)

def run(*args, **kwargs):
    return subprocess.check_call(*args, **kwargs)

def test_align_contigs():

    data_dir = join(data_base_dir, '0_RawData_PremadeGenome')
    exp_dir = join(data_base_dir, '1_alignContigs')
    out_dir = out_base_dir
    contig_dir = join(out_dir, 'premadeoutput')

    command = [
        'sisrs-python',
        '-p', '1',
        '-a', 'premade',
        '-c', '0',
        '-f', data_dir,
        '-z', out_dir,
        'align_contigs'
    ]
    run(command)

    # verify contigs.fa backed up
    assert(exists(join(contig_dir, 'contigs_OriginalNames.fa')))

    contig_file_path = join(contig_dir, 'contigs.fa')
    with open(contig_file_path, 'r') as f:
        first_line = f.next()

    # verify renaming
    assert(first_line.startswith('>SISRS_'))

    # index files generated
    index_filepaths = [
        join(contig_dir, filename) for filename in [
            'contigs.1.bt2',
            'contigs.2.bt2',
            'contigs.3.bt2',
            'contigs.4.bt2',
            'contigs.rev.1.bt2',
            'contigs.rev.2.bt2',
        ]
    ]

    for index_filepath in index_filepaths:
        assert(exists(index_filepath))

    


#def test_output_alignment():
#
#    data_dir = join(data_base_dir, '2_identifyFixedSites')
#    exp_dir = join(data_base_dir, '3_outputAlignment')
#    out_dir = out_base_dir
#
#    command = [
#        'sisrs-python',
#        '-f', data_dir,
#        '-z', out_dir,
#        'output_alignment'
#    ]
#    run(command)
#
#    assert cmp(
#        join(out_dir, 'alignment.nex'),
#        join(exp_dir, 'alignment.nex'))
