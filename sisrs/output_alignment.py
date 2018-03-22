from .command import Command
from .get_alignment import main as get_alignment


class OutputAlignmentCommand(Command):

    def run(self):
        args = self._args

        data_dir = args['data_dir']
        out_dir = args['out_dir']
        assembler = args['assembler']
        dir_lists = args['dir_lists']

        all_dirs = dir_lists.get_all_dirs()
        num_missing = len(all_dirs) - 2
        get_alignment(num_missing, 'X', out_dir, assembler)


