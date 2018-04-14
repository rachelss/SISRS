from abc import ABCMeta, abstractmethod


class Assembler(metaclass=ABCMeta):

    def __init__(self, dir_lists, out_dir, num_processors=1, **kwargs):
        self._num_processors = num_processors
        self._dir_lists = dir_lists
        self._out_dir = out_dir 

        self._kmer_size = kwargs['kmer_size']

    @abstractmethod
    def assemble():
        """Perform assembly"""
