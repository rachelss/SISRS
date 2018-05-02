import os
import subprocess
import shutil
import filecmp
import pickle
from filecmp import cmp
from os.path import join, exists

data_base_dir = 'pipeline_stages'
out_base_dir = 'output'
num_proc = 1

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

def setup():
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
    return result == b''

def pileups_match(f1, f2):

    return files_match(f1, f2)

def files_match(f1, f2):
    print("comparing {} and {}".format(f1, f2))
    return filecmp.cmp(f1, f2, shallow=False)

def pkls_match(path1, path2):
    print("comparing {} and {}".format(path1, path2))
    with open(path1, 'rb') as f1, open(path2, 'rb') as f2:
        p1 = pickle.load(f1)
        p2 = pickle.load(f2)
    return p1 == p2

def dirs_match(d1, d2):
    match, mismatch, errors = filecmp.cmpfiles(d1, d2, os.listdir(d1),
        shallow=False)

    return len(mismatch) == 0 and len(errors) == 0

def exists(f1):
    print("Checking if {} exists".format(f1))
    return os.path.exists(f1)


def test_align_contigs():

    setup()

    data_dir = join(data_base_dir, '0_RawData_PremadeGenome')
    exp_base_dir = join(data_base_dir, '1_alignContigs')
    contig_dirname = 'premadeoutput'
    contig_dir = join(out_base_dir, contig_dirname)

    command = [
        'sisrs-python',
        '-p', str(num_proc),
        '-a', 'premade',
        '-c', '0',
        '-f', data_dir,
        '-z', out_base_dir,
        'alignContigs'
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

    setup()

    data_dir = join(data_base_dir, '1_alignContigs')
    exp_base_dir = join(data_base_dir, '2_identifyFixedSites')
    contig_dirname = 'premadeoutput'
    contig_dir = join(out_base_dir, contig_dirname)

    command = [
        'sisrs-python',
        '-p', str(num_proc),
        '-a', 'premade',
        '-c', '0',
        '-f', data_dir,
        '-z', out_base_dir,
        'identifyFixedSites',
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

        assert(pkls_match(
            join(out_dir, 'pruned_dict.pkl'),
            join(exp_dir, 'pruned_dict.pkl')))

        filenames = [
            'contigs.1.bt2',
            'contigs.2.bt2',
            'contigs.3.bt2',
            'contigs.4.bt2',
            'contigs.fa',
            'contigs.fa.fai',
        ]
        for filename in filenames:
            assert(files_match(
                join(out_dir, filename),
                join(exp_dir, filename)))

        filenames = [
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

    setup()

    data_dir = join(data_base_dir, '2_identifyFixedSites')
    exp_dir = join(data_base_dir, '3_outputAlignment')
    out_dir = out_base_dir

    command = [
        'sisrs-python',
        '-f', data_dir,
        '-z', out_dir,
        '-c', '0',
        'outputAlignment'
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

    setup()

    data_dir = join(data_base_dir, '3_outputAlignment')
    exp_dir = join(data_base_dir, '4_changeMissing')
    out_dir = out_base_dir

    command = [
        'sisrs-python',
        '-f', data_dir,
        '-z', out_dir,
        'changeMissing',
        '-m', '0',
        '-c', '0',
    ]
    run(command)

    assert cmp(
        join(out_dir, 'alignment_m0.phylip-relaxed'),
        join(exp_dir, 'alignment_m0.phylip-relaxed'))

    assert cmp(
        join(out_dir, 'locs_m0.txt'),
        join(exp_dir, 'locs_m0.txt'))
    assert cmp(
        join(out_dir, 'locs_m0_Clean.txt'),
        join(exp_dir, 'locs_m0_Clean.txt'))


    assert exists(join(out_dir, 'alignmentDataWithoutSingletons'))
    assert cmp(
        join(out_dir, 'alignmentDataWithoutSingletons', 'alignment_pi.nex'),
        join(exp_dir, 'alignmentDataWithoutSingletons', 'alignment_pi.nex'))
    assert cmp(
        join(out_dir, 'alignmentDataWithoutSingletons', 'alignment_pi_m0.phylip-relaxed'),
        join(exp_dir, 'alignmentDataWithoutSingletons', 'alignment_pi_m0.phylip-relaxed'))
    assert cmp(
        join(out_dir, 'alignmentDataWithoutSingletons', 'locs_m0.txt'),
        join(exp_dir, 'alignmentDataWithoutSingletons', 'locs_m0.txt'))
    assert cmp(
        join(out_dir, 'alignmentDataWithoutSingletons', 'locs_m0_Clean.txt'),
        join(exp_dir, 'alignmentDataWithoutSingletons', 'locs_m0_Clean.txt'))


    assert exists(join(out_dir, 'alignmentDataWithOnlyBiallelic'))
    assert cmp(
        join(out_dir, 'alignmentDataWithOnlyBiallelic', 'alignment_bi.nex'),
        join(exp_dir, 'alignmentDataWithOnlyBiallelic', 'alignment_bi.nex'))
    assert cmp(
        join(out_dir, 'alignmentDataWithOnlyBiallelic', 'alignment_bi_m0.phylip-relaxed'),
        join(exp_dir, 'alignmentDataWithOnlyBiallelic', 'alignment_bi_m0.phylip-relaxed'))
    assert cmp(
        join(out_dir, 'alignmentDataWithOnlyBiallelic', 'locs_m0.txt'),
        join(exp_dir, 'alignmentDataWithOnlyBiallelic', 'locs_m0.txt'))
    assert cmp(
        join(out_dir, 'alignmentDataWithOnlyBiallelic', 'locs_m0_Clean.txt'),
        join(exp_dir, 'alignmentDataWithOnlyBiallelic', 'locs_m0_Clean.txt'))
