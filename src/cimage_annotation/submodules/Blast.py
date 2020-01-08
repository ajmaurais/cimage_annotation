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


def blastp(organism, database_path, query, verbose=False):

    # TODO: make blastp exe an option in setup.py
    exe = 'blastp'
    cmd = 'echo "{}"| {}'.format(query, exe)
    cmd += ' -query /dev/stdin'
    cmd += ' -db {}/{}'.format(database_path, DATABASES[organism])
    cmd += ' -outfmt ' + '0'
    cmd += ' -num_alignments ' + '5'
    cmd += ' -num_descriptions ' + '5'

    p = subprocess.Popen(cmd, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    if err and verbose:
        sys.stderr.write(err.decode('utf-8'))
    return p.returncode, out.decode('utf-8')

