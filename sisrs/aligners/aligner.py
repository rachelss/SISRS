from abc import ABCMeta, abstractmethod

class Aligner(metaclass=ABCMeta):

    def __init__(self, num_processors=1):
        self._num_processors = num_processors

    @abstractmethod
    def index():
        """Create index"""

    @abstractmethod
    def align():
        """Perform alignment"""
