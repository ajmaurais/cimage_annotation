
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

    assert(database_path is not None)
    _db_path = '{}/{}'.format(database_path, DATABASES[organism])

    # TODO: make blastp exe an option in setup.py
    exe = 'blastp'
    cmd = 'echo -e "{}"| {}'.format(query, exe)
    cmd += ' -query /dev/stdin'
    cmd += ' -db {}'.format(_db_path)
    cmd += ' -outfmt ' + '5'
    cmd += ' -num_alignments {}'.format(5)

    if verbose:
        sys.stdout.write('\n{}\n'.format(cmd))

    p = subprocess.Popen(cmd, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    if err and verbose:
        sys.stderr.write(err.decode('utf-8'))
    return p.returncode, out.decode('utf-8')

