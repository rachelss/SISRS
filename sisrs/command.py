from abc import ABCMeta, abstractmethod

class Command(metaclass=ABCMeta):

    def __init__(self, args):
        self._args = args

    @abstractmethod
    def run():
        """Run command"""
