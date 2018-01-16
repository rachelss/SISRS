from subprocess import Popen


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


