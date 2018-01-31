import os
import subprocess
import shutil
import filecmp
from filecmp import cmp
from os.path import join, exists

data_base_dir = 'pipeline_stages'
out_base_dir = 'output'
num_proc = '1'

taxon_names = [
    'GorGor',
    'HomSap',
    'HylMol',
    'MacFas',
    'MacMul',
    'PanPan',
    'PanTro',
    'PonPyg'
]

if os.path.exists(out_base_dir):
    shutil.rmtree(out_base_dir)

def run(*args, **kwargs):
    return subprocess.check_call(*args, **kwargs)

def bam_match(f1, f2):

    bamdiff_command = [
        'bamdiff',
        f1,
        f2
    ]
    
    print("comparing {} and {}".format(f1, f2))
    result = subprocess.check_output(bamdiff_command)
    return result == ''

def bam_match_old(f1, f2):
    bamhash_command_f1 = [
        'bamhash_checksum_bam',
        '--no-paired',
        f1,
    ]

    bamhash_command_f2 = [
        'bamhash_checksum_bam',
        '--no-paired',
        f2,
    ]
    
    output_f1 = subprocess.check_output(bamhash_command_f1)
    output_f2 = subprocess.check_output(bamhash_command_f2)

    print("comparing {} and {}".format(f1, f2))
    print(output_f1)
    print(output_f2)

    return output_f1 == output_f2

def pileups_match(f1, f2):

    return files_match(f1, f2)

def files_match(f1, f2):
    print("comparing {} and {}".format(f1, f2))
    return filecmp.cmp(f1, f2, shallow=False)

def dirs_match(d1, d2):
    match, mismatch, errors = filecmp.cmpfiles(d1, d2, os.listdir(d1),
        shallow=False)

    return len(mismatch) == 0 and len(errors) == 0

def exists(f1):
    print("Checking if {} exists".format(f1))
    return os.path.exists(f1)

def test_align_contigs():

    data_dir = join(data_base_dir, '0_RawData_PremadeGenome')
    exp_base_dir = join(data_base_dir, '1_alignContigs')
    contig_dirname = 'premadeoutput'
    contig_dir = join(out_base_dir, contig_dirname)

    command = [
        'sisrs-python',
        '-p', num_proc,
        '-a', 'premade',
        '-c', '0',
        '-f', data_dir,
        '-z', out_base_dir,
        'align_contigs'
    ]
    run(command)

    assert(dirs_match(
        join(out_base_dir, 'premadeoutput'),
        join(exp_base_dir, 'premadeoutput')))

    for taxon_name in taxon_names:
        out_dir = join(out_base_dir, taxon_name)
        exp_dir = join(exp_base_dir, taxon_name)

        out_bam_path = join(out_dir, taxon_name + '.bam')
        exp_bam_path = join(exp_dir, taxon_name + '.bam')

        assert(bam_match(out_bam_path, exp_bam_path))
        

def test_identify_fixed_sites():

    data_dir = join(data_base_dir, '1_alignContigs')
    exp_base_dir = join(data_base_dir, '2_identifyFixedSites')
    contig_dirname = 'premadeoutput'
    contig_dir = join(out_base_dir, contig_dirname)

    command = [
        'sisrs-python',
        '-p', num_proc,
        '-a', 'premade',
        '-c', '0',
        '-f', data_dir,
        '-z', out_base_dir,
        'identify_fixed_sites',
    ]
    run(command)

    assert(dirs_match(
        join(out_base_dir, 'premadeoutput'),
        join(exp_base_dir, 'premadeoutput')))

    for taxon_name in taxon_names:
        out_dir = join(out_base_dir, taxon_name)
        exp_dir = join(exp_base_dir, taxon_name)

        out_bam_path = join(out_dir, taxon_name + '.bam')
        exp_bam_path = join(exp_dir, taxon_name + '.bam')

        assert(bam_match(out_bam_path, exp_bam_path))


        assert(files_match(
            join(out_dir, 'pruned_dict.pkl'),
            join(exp_dir, 'pruned_dict.pkl')))

        filenames = [
            'contigs.1.bt2',
            'contigs.2.bt2',
            'contigs.3.bt2',
            'contigs.4.bt2',
            'contigs.fa',
            'contigs.fa.fai',
            taxon_name + '.pileups',
            taxon_name + '.bam.bai'
        ]
        for filename in filenames:
            assert exists(join(out_dir, filename))

    # TODO: see if we can figure out a way to compare pileup files
    #for taxon_name in taxon_names:
    #    out_dir = join(out_base_dir, taxon_name)
    #    exp_dir = join(exp_base_dir, taxon_name)

    #    out_pileup_path = join(out_dir, taxon_name + '.pileups')
    #    exp_pileup_path = join(exp_dir, taxon_name + '.pileups')

    #    assert(pileups_match(out_pileup_path, exp_pileup_path))


def test_output_alignment():

    data_dir = join(data_base_dir, '2_identifyFixedSites')
    exp_dir = join(data_base_dir, '3_outputAlignment')
    out_dir = out_base_dir

    command = [
        'sisrs-python',
        '-f', data_dir,
        '-z', out_dir,
        'output_alignment'
    ]
    run(command)

    assert cmp(
        join(out_dir, 'alignment.nex'),
        join(exp_dir, 'alignment.nex'))

    assert cmp(
        join(out_dir, 'alignment_bi.nex'),
        join(exp_dir, 'alignment_bi.nex'))

    assert cmp(
        join(out_dir, 'alignment_pi.nex'),
        join(exp_dir, 'alignment_pi.nex'))

def test_change_missing():

    data_dir = join(data_base_dir, '3_outputAlignment')
    exp_dir = join(data_base_dir, '4_changeMissing')
    out_dir = out_base_dir

    command = [
        'sisrs-python',
        '-f', data_dir,
        '-z', out_dir,
        'change_missing',
    ]
    run(command)

    assert cmp(
        join(out_dir, 'alignment_m6.phylip-relaxed'),
        join(exp_dir, 'alignment_m6.phylip-relaxed'))

    assert cmp(
        join(out_dir, 'locs_m6.txt'),
        join(exp_dir, 'locs_m6.txt'))

    assert cmp(
        join(out_dir, 'locs_m6_Clean.txt'),
        join(exp_dir, 'locs_m6_Clean.txt'))
