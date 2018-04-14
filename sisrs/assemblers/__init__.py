from .velvet import VelvetAssembler
import sys

def create_assembler(dir_lists, out_dir, assembler_type='velvet', **kwargs):

    if assembler_type == 'velvet':
        return VelvetAssembler(dir_lists, out_dir, **kwargs)
