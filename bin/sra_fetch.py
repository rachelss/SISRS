#!/usr/bin/env python
"""
Fetch all FASTQ files associated with an SRA ID with a SISRS-friendly directory structure.

usage:

    sra_fetch.py [sra_id] [url]

"""
import urllib
import urllib2
import os
import shutil
import sys
import re
import time
from ftplib import FTP

CHUNK = 16 * 1024
file_suffix = "{0}.fastq.gz"
base_dl_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/{1}/{2}.fastq.gz"
base_species_url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=read_run&fields=scientific_name"

#Gets the species for the given SRA_ID
def fetch_species(SRA_ID):
    file = urllib2.urlopen(base_species_url.format(SRA_ID))
    return file.readlines()[1].strip()

def download_file(filename, filedest):
    response = urllib2.urlopen(filename)
    with open(filedest, "wb") as out:
        while True:
            next_chunk = response.read(CHUNK)
            if not next_chunk:
                response.close()
                time.sleep(1)#let's not DoS the ENA
                break
            out.write(next_chunk)

def main():
    SRA_ID = sys.argv[1]
    URL_LIST = sys.argv[2]
    directory = re.sub(r' ',r'_',fetch_species(SRA_ID))
    directory = re.sub(r'-',r'_',directory)
    directory = "sra_fastqs/" + directory
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = directory + "/"
    sra_url = base_dl_url.format(SRA_ID[:6],SRA_ID,SRA_ID)
    s = URL_LIST.split(';',1)

    if len(s) == 2:
        s[1]=s[1].strip(';')
        i = 0;
        for u in s:
            i = i + 1
            filename = "ftp://" + u
            filedest = directory + file_suffix.format(SRA_ID + "_R" + str(i))
            download_file(filename, filedest)
    else:
        filename = "ftp://" + URL_LIST
        filedest = directory + file_suffix.format(SRA_ID)
        download_file(filename, filedest)


if __name__ == "__main__":
    main()
