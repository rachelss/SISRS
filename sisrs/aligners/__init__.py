from .bowtie2 import Bowtie2Aligner
import sys

def create_aligner(aligner_type='bowtie2', **kwargs):

    if aligner_type == 'bowtie2':
        return Bowtie2Aligner(**kwargs)
