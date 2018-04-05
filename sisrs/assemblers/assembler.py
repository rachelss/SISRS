from abc import ABCMeta, abstractmethod


class Assembler(metaclass=ABCMeta):

    def __init__(self, num_processors=1):
        self._num_processors = num_processors

    @abstractmethod
    def assemble():
        """Perform assembly"""
