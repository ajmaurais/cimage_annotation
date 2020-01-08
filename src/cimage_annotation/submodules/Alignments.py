#! usr/bin/env/ python

import sys
import re
from multiprocessing import Pool
from multiprocessing import cpu_count
import itertools
from tqdm import tqdm
import functools

from .Blast import blastp

UNIPROT_RE = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'

def align_write(fname, dat):
    with open(fname, 'a') as outF:
        outF.write(dat + '\n\n\n\n\n\n')


def align_file_parse(position, path=None, string=None):

    if sum([bool(x) for x in [path, string]]) != 1:
        raise RuntimeError('path xor string must be specified!')

    # get arg for text to process
    if path is not None:
        f = open(path + 'alignment.txt', 'r')
        align = f.read()
        f.close()
    elif string is not None:
        align=string

    conserved = ''
    id = ''
    evalue = ''
    homolog_position = ''
    seq_list = []
    if len(align.split('\n\n')) < 9:
        conserved = 'Error'
    else:
        align_header = align.split('\n\n')[7]

        if align_header.strip() == '***** No hits found *****':
            conserved = '--'
        else:
            aligned_proteins = align_header.split('\n')

            r = '^([ts][rp]\|){0,1}(' + UNIPROT_RE + ').*\s+([\d\.\-e]+)\s*$'
            match = re.search(r, aligned_proteins[0])
            if not match:
                raise RuntimeError('Failed to parse aligned_proteins!\n\tFailed on value: {}'.format(aligned_proteins[0]))
            id = match.group(2)
            evalue = match.group(4)

            if float(evalue) > float('1e-5'):
                conserved = '--'
            else:

                align_data = align.split('>')[1].split('\n\n\n')[0].split('\n\n')

                for item in align_data:
                    if item[:5] == 'Query':
                        query = item.split('\n')[0].split()
                        subject = item.split('\n')[2].split()

                        if int(position) >= int(query[1]) and int(position) <= int(query[3]):
                            m = int(query[1])
                            n = int(subject[1])
                            i = 0
                            while i < len(query[2]):
                                if query[2][i] == '-':
                                    a = m
                                else:
                                    a = m
                                    m += 1
                                if subject[2][i] == '-':
                                    b = n
                                else:
                                    b = n
                                    n += 1
                                seq_list.append(str(a) + '\t' + query[2][i] + '\t' + str(b) + '\t' + subject[2][i])

                                i += 1
                        else:
                            pass
                    else:
                        pass

            if seq_list == []:
                conserved = '--'
            else:
                for item in seq_list:
                    if item.split('\t')[0] == str(position) and item.split('\t')[1] == item.split('\t')[3]:
                        conserved = 'Yes'
                        homolog_position = item.split('\t')[2]
                    else:
                        pass
            if homolog_position == '' and conserved == '':
                conserved = 'No'

    return id, evalue, conserved, homolog_position


def _blastp_worker(search_item, db = None, verbose=False):

    return_code, dat = blastp(search_item[2], db, search_item[1], verbose=verbose)
    return dat


def align_all(peptides, sequences, db_path, organisms, nThread=None, show_bar=True, verbose=False):

    #construct list to pass to blastp worker
    search_list = list()
    for id in set([x['id'] for x in peptides]):
        for o in organisms:
            search_list.append((id, sequences[id], o))

    #calculate number of threads required
    _nThread = int(1)
    listLen = len(search_list)
    cpuCount = cpu_count()
    if nThread is None:
        _nThread = cpuCount if cpuCount < listLen else listLen
    else:
        _nThread = nThread

    sys.stdout.write('Performing alignment with {} thread(s)...\n'.format(_nThread))
    results = list()
    if show_bar:
        with Pool(processes=_nThread) as pool:
            results = list(tqdm(pool.imap(functools.partial(_blastp_worker, db = db_path, verbose=verbose),
                                          search_list),
                                          total = listLen,
                                          miniters=1,
                                          file = sys.stdout))
    else:
        for it in search_list:
            results.append(_blastp_worker(it, db=db_path, verbose=verbose))

    assert(len(search_list) == len(results))

    ret=dict()
    for sl, r in zip(search_list, results):
        if sl[0] not in ret:
            ret[sl[0]]=dict()
        ret[sl[0]][sl[2]]=r

    return ret

