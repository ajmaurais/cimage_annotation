#! /usr/bin/env/ python

import sys
import os
import os.path
from urllib.error import URLError, HTTPError
from http.client import BadStatusLine
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import cpu_count
from tqdm import tqdm
import functools
from Bio import ExPASy, SeqIO, SwissProt


features_list =['CA_BIND', 'ZN_FING', 'DNA_BIND', 'NP_BIND',
                 'ACT_SITE', 'METAL', 'BINDING', 'SITE',
                 'NON_STD', 'MOD_RES', 'LIPID', 'CARBOHYD',
                 'DISULFID', 'CROSSLINK', 'VARIANT', 'MUTAGEN',
                 'UNSURE', 'CONFLICT', 'REGION']


def _parse_record(handle):

    record = None
    if handle is not None:
        try:
            record = SwissProt.read(handle)
        except ValueError as e:
            pass
    return record


def make_request(uniprot_id, verbose = True, n_retry = 10):
    '''
    ExPASy get_sprot_raw wrapper to make retries if an http error occurs.

    Paramaters
    ----------
    uniprot_id: str
        Uniprot ID to retreive.

    n_retry: int
        Number of times to retry request if an error occurs
    '''

    n_iter = n_retry if n_retry > 0 else 1
    ret = None
    for i in range(n_iter):
        try:
            handle = ExPASy.get_sprot_raw(uniprot_id)
        except HTTPError as e:
            if e.getcode() >= 400:
                if verbose:
                    sys.stderr.write('No UniProt page found for {}\n\tStatus code: {}\n'.format(uniprot_id, e.getcode()))
                break
            else:
                if verbose:
                    sys.stderr.write('Retry {} of {} for {}\n\t{}\n'.format(i, n_iter, uniprot_id, e))
                continue
        except (BadStatusLine, URLError) as e:
            if verbose:
                sys.stderr.write('Retry {} of {} for {}\n\t{}\n'.format(i, n_iter, uniprot_id, e))
            continue
        else:
            return _parse_record(handle)

    return ret


def get_uniprot_records(ids, parallel, verbose=False):
    '''
    Get a dict of UniProt records.

    Parameters
    ----------
    ids: list like
        List of uniprot ids to retreive.
    parallel: bool
        Should the http request be made in parallel?

    Return
    ------
    records: dict
        Key value paris of IDs and records.
    '''

    #calculate number of threads required
    if parallel:
        _nThread = int(1)
        listLen = len(ids)
        cpuCount = cpu_count()
        _nThread = cpuCount if cpuCount < listLen else listLen
    else:
        _nThread = 1

    sys.stdout.write('Searching for data with {} thread(s)...\n'.format(_nThread))
    ret = list()
    with Pool(processes=_nThread) as pool:
        ret = list(tqdm(pool.imap(functools.partial(make_request, verbose=verbose), ids),
                             total = listLen,
                             miniters=1,
                             file = sys.stdout))

    assert(len(ids) == len(ret))
    return {k: record for k, record in zip(ids, ret)}

def protein_location(record):
    pro_location = ''
    for item in record.comments:

        if item.split(':')[0] == 'SUBCELLULAR LOCATION':
            if item.split(':')[1].strip().split(' ')[0] == 'Isoform':
                pro_location += item.split(':')[2]
            else:
                pro_location += item.split(':')[1]
        else:
            pass
    return str(pro_location)


def cys_position(record, sequence, location):
    peptide_location = str(record.sequence).find(sequence)
    if str(peptide_location) == str(-1):
        cys_location = 0
    else:
        cys_location = peptide_location + location

    return str(cys_location)


def cys_function(record, position):
    cys_function = ''

    for feature in record.features:
        try:
            if feature.type.upper() == 'DISULFID' and \
               str(feature.location.start)[0] != '?' and str(feature.location.end)[0] != '?' and \
               int(position) >= int(feature.location.start) and int(position) <= int(feature.location.end):
                cys_function += (str(feature.type) + '--' + str(feature.qualifiers) + ' || ')

            elif feature.type.upper() in features_list and \
                 str(feature.location.start)[0] != '?' and str(feature.location.end)[0] != '?' and \
                 (int(feature.location.end) - int(feature.location.start)) <= 10 and \
                 int(position) >= int(feature.location.start) and int(position) <= int(feature.location.end):
                cys_function += (str(feature.type) + '--' + str(feature.qualifiers) + ' || ')
        except TypeError as e:
            continue

    return cys_function

def ExPasy(id, sequence, location, record):
    position = ''
    function = ''
    organism = ''
    full_sequence = ''
    pro_location = ''

    if record is not None:
        organism = record.organism
        pro_location = protein_location(record)
        position = cys_position(record, sequence, location)
        function = cys_function(record, position)
        full_sequence = record.sequence
    else:
        position = 'Bad ID'
    return organism, position, function, full_sequence, pro_location

def ExPasy_alt(id, position, record):
    protein = ''
    function = ''

    if record is not None:
        protein = record.description
        function = cys_function(record, position)

    else:
        protein = 'Bad ID'

    return protein, function

