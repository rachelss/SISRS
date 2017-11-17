import os
import subprocess
import shutil
from filecmp import cmp
from os.path import join

data_base_dir = 'test_data'
out_base_dir = 'output'

print("out:")
print(out_base_dir)
if os.path.exists(out_base_dir):
    print("remove it")
    shutil.rmtree(out_base_dir)

def run(*args, **kwargs):
    return subprocess.check_call(*args, **kwargs)


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
