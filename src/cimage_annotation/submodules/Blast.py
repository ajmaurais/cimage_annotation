#! /usr/bin/env/ python

import sys
import os
import subprocess


DATABASES = {'human':'human_nr_uniprot',
             'mouse':'mouse_nr_uniprot',
             'fly':'fly_nr_uniprot',
             'yeast':'yeast_nr_uniprot',
             'mustard':'mustard_nr_uniprot',
             'worms':'worms_nr_uniprot'}


def blastp(organism, path, database_path):

    exe = 'blastp'
    cmd = exe
    cmd += ' -query ' + path + 'sequence.txt'
    cmd += ' -db {}/{}'.format(database_path, DATABASES[organism])
    cmd += ' -out ' + path + 'alignment.txt'
    cmd += ' -outfmt ' + '0'
    cmd += ' -num_alignments ' + '5'
    cmd += ' -num_descriptions ' + '5'

    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)
    return sts[1] # exit status
