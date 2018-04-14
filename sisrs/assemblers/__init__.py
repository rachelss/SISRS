from .velvet import VelvetAssembler
from .minia import MiniaAssembler
import sys

def create_assembler(assembler_type, dir_lists, out_dir, **kwargs):

    if assembler_type == 'velvet':
        return VelvetAssembler(dir_lists, out_dir, **kwargs)
    elif assembler_type == 'minia':
        return MiniaAssembler(dir_lists, out_dir, **kwargs)
    else:
        raise Exception("Invalid assembler type: {}".format(assembler_type))
