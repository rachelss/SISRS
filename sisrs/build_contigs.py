from .command import Command
from .assemblers import create_assembler


class BuildContigsCommand(Command):
    
    def run(self):

        assembler = create_assembler(
                self._args['assembler'],
                self._args['dir_lists'],
                self._args['out_dir'],
                num_processors=self._args['num_processors'],
                ref_file_path=self._args['reference'],
                kmer_size=self._args['kmer_size'])

        assembler.assemble()
