import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
    name='sisrs',
    version='0.1.0',
    packages=find_packages(),
    description="Site Identification from Short Read Sequences.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://cartwrig.ht/software/',
    license='GPLv3',
    install_requires=[
        'biopython',
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'sisrs-python = sisrs:main'
        ],
    },
    keywords='bioinformatics genetics',
)
