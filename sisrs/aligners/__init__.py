from .bowtie2 import Bowtie2Aligner

def create_aligner(aligner_type='bowtie2'):

        if aligner_type == 'bowtie2':
            return Bowtie2Aligner()
