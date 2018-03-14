import os
from subprocess import Popen
from glob import glob
from subprocess import Popen, PIPE 


class Process(object):

    def __init__(self, command, stdin=None, stdout=None):
        try:
            self.proc = Popen(command, stdin=stdin, stdout=stdout)
        except Exception as e:
            print("Error calling: {}".format(command[0]))
            raise e 

    def wait(self):
        return self.proc.wait()

    def pipe(self):
        return self.proc.stdout


class MpileupProcess(object):

    def __init__(self, dir_, contig_file_path, num_processors):

        taxon_name = os.path.basename(dir_)
        bam_path = os.path.join(dir_, '{}.bam'.format(taxon_name))
        pileups_path  = os.path.join(dir_, '{}.pileups'.format(taxon_name))

        command = [
            'samtools',
            '-f', contig_file_path,
            bam_path,
            '-o', pileups_path
        ]
        Process(command).wait()

    #parallel --jobs "${PROCESSORS}" 'samtools mpileup -f' "${OUTFOLDER}"/"${CONTIGS}"/contigs.fa '"$( echo {}/$(basename {} ) )".bam' '> "$( echo {}/$(basename {} ) )".pileups' ::: "${FOLDERLISTA[@]}"
