from .velvet import VelvetAssembler
import sys

def create_assembler(assembler_type='velvet', **kwargs):

    if assembler_type == 'velvet':
        return VelvetAssembler(**kwargs)
