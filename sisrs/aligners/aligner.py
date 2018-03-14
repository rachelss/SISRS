from abc import ABCMeta, abstractmethod

class Aligner(metaclass=ABCMeta):

    def __init__(self):
        """Base constructor"""

    @abstractmethod
    def index():
        """Create index"""

    @abstractmethod
    def align():
        """Perform alignment"""
