import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
    # TODO: upgrade to python 3
#with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
    name='sisrs',
    version='0.0.1',
    packages=find_packages(),
    description="Site Identification from Short Read Sequences.",
    long_description=long_description,
    url='https://cartwrig.ht/software/',
    install_requires=[
        'biopython',
        'click',
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'sisrs-python = sisrs:main'
        ],
    },
)
