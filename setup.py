import os
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='sisrs',
    version='0.0.1',
    description="Site Identification from Short Read Sequences.",
    long_description=long_description,
    url='https://cartwrig.ht/software/',
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'sisrs-python = sisrs:main'
        ],
    },
    scripts=['bin/sisrs']
)
