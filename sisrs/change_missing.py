import os

from filter_nexus_for_missing import main as filter_nexus_for_missing
from process import Process 
from subprocess import PIPE 

def clean_file(in_path, out_path):

    grep_command = [
        'grep', '-oe', "SISRS_[^/]*",
        in_path,
    ]
    grep_proc = Process(grep_command, stdout=PIPE)

    uniq_command = [
        'uniq', '-c',
    ]
    uniq_proc = Process(uniq_command, stdin=grep_proc.pipe(), stdout=PIPE)

    sort_command = [
        'sort', '-k1', '-nr'
    ]
    sort_proc = Process(sort_command, stdin=uniq_proc.pipe(), stdout=PIPE)

    awk_command = [
        'awk', "{print $2}",
    ]
    awk_proc = Process(awk_command, stdin=sort_proc.pipe(), stdout=PIPE)

    # tee is simply for copying to file, since apparently awk can only deal
    # with stdout
    tee_command = [
        'tee', out_path
    ]
    tee_proc = Process(tee_command, stdin=awk_proc.pipe(), stdout=PIPE)

    tee_proc.wait()


    # TODO: fix python implementation
    # For some reason the python code below has some subtle differences in
    # the output files at the binary level. The text is the same but the hex
    # dumps are different

    #with open(in_path, 'r') as f_in:
    #    
    #    lines = [ line[:line.find('/')] for line in f_in ]

    #counts = {}
    #for name in lines:
    #    if name not in counts:
    #        counts[name] = 0
    #    counts[name] += 1

    #print(counts)

    #sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    #with open(out_path, 'w') as f_out:

    #    for count in sorted_counts:
    #        f_out.write("{}\n".format(count[0]))


class ChangeMissingCommand(object):

    def __init__(self, data, missing):

        self._data = data 
        self._missing = missing

    def run(self):

        data = self._data
        missing = self._missing

        alignment_file_path = os.path.join(data['out_dir'], 'alignment.nex')


        if missing is None:
            dir_lists = data['dir_lists']
            all_dirs = dir_lists.get_all_dirs()
            str_missing = str(len(all_dirs) - 2)
        else:
            str_missing = str(missing)

        filter_nexus_for_missing(alignment_file_path, str_missing)

        locs_filename = 'locs_m{}.txt'.format(str_missing)
        locs_file_path = os.path.join(data['out_dir'], locs_filename)

        locs_clean_filename = 'locs_m{}_Clean.txt'.format(str_missing)
        locs_clean_file_path = os.path.join(
            data['out_dir'], locs_clean_filename)

        clean_file(locs_file_path, locs_clean_file_path)
