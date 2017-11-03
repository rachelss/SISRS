import subprocess
import sys

def main():
    command = ['sisrs'] +  sys.argv[1:]
    subprocess.run(command)
