
import sys
import re
from urllib.error import URLError
from http.client import BadStatusLine
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import cpu_count
from tqdm import tqdm
import functools
from Bio import ExPASy, SwissProt


features_list = ['CA_BIND', 'ZN_FING', 'DNA_BIND', 'NP_BIND',
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


def make_request(uniprot_id, verbose=True, n_retry=10):
    '''
    ExPASy get_sprot_raw wrapper to make retries if an http error occurs.

    Parameters
    ----------
    uniprot_id: str
        Uniprot ID to retrieve.

    n_retry: int
        Number of times to retry request if an error occurs
    '''

    n_iter = n_retry if n_retry > 0 else 1
    ret = None
    for i in range(n_iter):
        try:
            handle = ExPASy.get_sprot_raw(uniprot_id)
        except ValueError as e:
            if verbose:
                sys.stderr.write('No UniProt page found for {}\n'.format(uniprot_id))
            break
        except (BadStatusLine, URLError) as e:
            if verbose:
                sys.stderr.write('Retry {} of {} for {}\n\t{}\n'.format(i, n_iter, uniprot_id, e))
            continue
        else:
            return _parse_record(handle)

    return ret


def get_uniprot_records(ids, nThread, verbose=False, show_bar=True):
    '''
    Get a dict of UniProt records.

    Parameters
    ----------
    ids: list like
        List of uniprot ids to retrieve.
    nThread: int
        Specify how many threads to use to retrieve Uniprot records.
    verbose: bool
        Verbose output?
    show_bar: bool
        Shouold status bar be shown?

    Return
    ------
    records: dict
        Key value paris of IDs and records.
    '''

    #calculate number of threads required
    _nThread = int(1)
    listLen = len(ids)
    cpuCount = cpu_count()
    if nThread is None:
        _nThread = cpuCount if cpuCount < listLen else listLen
    else:
        _nThread = nThread

    sys.stdout.write('Searching for data with {} thread(s)...\n'.format(_nThread))
    ret = list()
    if show_bar:
        with Pool(processes=_nThread) as pool:
            ret = list(tqdm(pool.imap(functools.partial(make_request, verbose=verbose), ids),
                                 total = listLen,
                                 miniters=1,
                                 file = sys.stdout))
    else:
        length = len(ids)
        for i, it in enumerate(ids):
            sys.stdout.write('Working on {} of {}\n'.format(i, length))
            ret.append(make_request(it, verbose=verbose))

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
    return str(pro_location).strip()


def cys_position(protein_seq, peptide_seq, mod_loc):
    '''
    Get the position of a cysteine in a parent protein.

    Parameters
    ----------
    protein_seq: str
        Full parent protein sequence.
    peptide_seq: str
        Peptide sequence (Modifications are removed automatically if they are present.)
    mod_loc: int
        Index of cysteine position in parent protein.

    Returns
    -------
    cys_loc: int
        Cysteine location in parent protein. (Starting from index 0.)
    '''

    _peptide_seq = peptide_seq.replace('*', '')
    peptide_start = str(protein_seq).find(_peptide_seq)
    if str(peptide_start) == str(-1):
        cys_loc = -1
    else:
        cys_loc = peptide_start + mod_loc

    return cys_loc


def res_features(record, position, all_features=False):
    '''
    Get residue function annotation at `position` if it exists in `record`.

    Parameters
    ----------
    position: int
        Index (starting from 0) of the residue of interest.
    record: Bio.SwissProt.Record
        Record to get data from.
    all_features: bool
        Should all features or simplified list be used?

    Returns
    -------
    res_features: str
        Annotaton for cysteine function.
    '''

    ret = ''
    domains = ''
    for feature in record.features:
        try:
            if feature.type.upper() == 'DOMAIN' and \
                 str(feature.location.start)[0] != '?' and str(feature.location.end)[0] != '?' and \
                 int(position) >= int(feature.location.start) and int(position) < int(feature.location.end):
                     domains += (str(feature.type) + '--' + str(feature.qualifiers) + ' || ')
            if all_features:
                if str(feature.location.start)[0] != '?' and str(feature.location.end)[0] != '?' and \
                   int(position) >= int(feature.location.start) and int(position) < int(feature.location.end):
                    ret += (str(feature.type) + '--' + str(feature.qualifiers) + ' || ')
            else:
                if feature.type.upper() == 'DISULFID' and \
                   str(feature.location.start)[0] != '?' and str(feature.location.end)[0] != '?' and \
                   int(position) >= int(feature.location.start) and int(position) < int(feature.location.end):
                    ret += (str(feature.type) + '--' + str(feature.qualifiers) + ' || ')

                elif feature.type.upper() in features_list and \
                     str(feature.location.start)[0] != '?' and str(feature.location.end)[0] != '?' and \
                     (int(feature.location.end) - int(feature.location.start)) <= 10 and \
                     int(position) >= int(feature.location.start) and int(position) < int(feature.location.end):
                    ret += (str(feature.type) + '--' + str(feature.qualifiers) + ' || ')
        except TypeError as e:
            continue

    return ret, domains

def parse_domains(domain_s):
    ret = ''
    m = re.findall(r'\'note\': \'([\w\s\-_]+)\'', domain_s)
    if m:
        ret = '|'.join(set(m))
    return ret

def ExPasy(sequence, record, all_features=False, res_sep='|', fxn_sep='!', combine_method=1):

    if combine_method != 1:
        raise NotImplementedError('combine_method {} not implemented.'.format(combine_method))

    seq_no_mod = sequence.replace('*', '')
    position = ''
    function = ''
    domain = ''
    organism = ''
    full_sequence = ''
    pro_location = ''

    if record is not None:
        organism = record.organism
        pro_location = protein_location(record)
        full_sequence = record.sequence

        positions = list()
        functions = list()
        domains = list()
        for i, x in enumerate(re.finditer('\*', sequence)):
            mod_loc = x.start()-(i+1)
            cys_pos = cys_position(full_sequence, seq_no_mod, mod_loc)
            if cys_pos == -1:
                positions.append('RESIDUE_NOT_FOUND')
            else:
                positions.append(str(cys_pos + 1)) # convert to 1 based indexing here
                functions_t, domains_t = res_features(record, cys_pos, all_features=all_features)
                functions.append(functions_t)
                domains.append(domains_t)

        position = res_sep.join(positions)
        if ''.join(functions):
            function = functions[0] if len(functions) == 1 else fxn_sep.join(['{}:{}'.format(p, s) for p, s in zip(positions, functions)])
        if ''.join(domains):
            domain = domains[0] if len(domains) == 1 else fxn_sep.join(['{}:{}'.format(p, s) for p, s in zip(positions, domains)])
            domain = parse_domains(domain)

    else:
        position = 'BAD_ID'
    return organism, position, function, domain, full_sequence, pro_location

