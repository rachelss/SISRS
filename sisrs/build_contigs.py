from .command import Command
from .assemblers import create_assembler


class BuildContigsCommand(Command):
    
    def run(self):
        print("run BuildContigsCommand")
        print(self._args['assembler'])

        assembler = create_assembler(self._args['assembler'])
        assembler.assemble()
