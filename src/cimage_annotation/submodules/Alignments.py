
import sys
import re
import xml.etree.ElementTree as ET
from multiprocessing import Pool
from multiprocessing import cpu_count
import itertools
from tqdm import tqdm
import functools

from .Blast import blastp


class Alignment(object):
    '''
    Read BLAST XML file and provide methods to acess underlying alignment data.
    '''

    def __init__(self, raw_xml, query_id = None, query_description = None, query_organism = None):
        '''
        Default constructor

        Parameters
        ----------
        raw_xml: str
            BLAST XML file as text.
        '''

        self.query_id = query_id
        self.query_description = query_description
        self.query_organism = query_organism

        if raw_xml:
            self._tree = ET.fromstring(raw_xml)
            self._best_hit = self._tree.find('.//Iteration_hits/Hit')
            self._hsp = None
            if self._best_hit is None:
                self._empty = True
            else:
                self._empty = False
        else:
            self._empty = True


    def _get_hit_path_text(self, path):
        temp = self._best_hit.find(path)
        return None if temp is None else temp.text


    def get_best_id(self):
        '''
        Get acession of best alignmant match.

        Returns
        -------
        best_id: str
        '''
        if self._empty:
            return ''
        return self._get_hit_path_text('./Hit_accession')


    def get_best_description(self):
        '''
        Get description of best alignmant match.

        Returns
        -------
        best_description: str
        '''
        if self._empty:
            return ''
        return self._get_hit_path_text('./Hit_def')


    def get_best_evalue(self):
        '''
        Get evalue of best alignment.

        Returns
        -------
        evalue: float
        '''
        if self._empty:
            return None
        return float(self._get_hit_path_text('./Hit_hsps/Hsp/Hsp_evalue'))


    def _populate_hsp(self):
        self._hsp = {x.tag: x.text for x in self._best_hit.findall('./Hit_hsps/Hsp/')}

        self._hsp['hit_seq_map'] = list()
        seq_index = 0
        for i, c in enumerate(self._hsp['Hsp_hseq']):
            self._hsp['hit_seq_map'].append(seq_index)
            if re.match('[A-Za-z]', c):
                seq_index += 1

        self._hsp['query_seq_map'] = list()
        for i, c in enumerate(self._hsp['Hsp_qseq']):
            if re.match('[A-Za-z]', c):
                self._hsp['query_seq_map'].append(i)

        # make the things which need to be numbers, numbers
        for entry in ('Hsp_hit-to', 'Hsp_hit-from', 'Hsp_query-to', 'Hsp_query-from'):
            self._hsp[entry] = int(self._hsp[entry])
        self._hsp['Hsp_evalue'] = float(self._hsp['Hsp_evalue'])


    def conserved_at_position(self, pos):
        '''
        Determine whether residue at `pos` is conserved in hit sequence.

        Returns
        -------
        conserved: bool
            True if conserved.
        '''

        if self._empty:
            return False

        # populate self._hsp with dict of dat from best hit
        if self._hsp is None:
            self._populate_hsp()

        if pos < self._hsp['Hsp_query-from'] or pos > self._hsp['Hsp_query-to']:
            return False

        align_index = pos - self._hsp['Hsp_query-from']
        query_index = self._hsp['query_seq_map'][align_index]
        if query_index == -1:
            return False
        query_residue = self._hsp['Hsp_qseq'][query_index]
        hit_residue = self._hsp['Hsp_hseq'][query_index]
        sys.stdout.write('{}: {} == {} -> {}\n'.format(pos, query_residue, hit_residue, (hit_residue == query_residue)))

        return hit_residue == query_residue


    def alignment_at_position(self, pos):
        '''
        Get alignment at specified position in query sequence.

        Parameters
        ----------
        pos: int
            Index of query position of interest. (starting from 1)

        Returns
        -------
        residue, position: str, int
            Residue at `pos` and position in aligned sequence.
            If no match, returns None, None.
        '''

        if self._empty:
            return None, None

        # populate self._hsp with dict of dat from best hit
        if self._hsp is None:
            self._populate_hsp()

        if pos < self._hsp['Hsp_query-from'] or pos > self._hsp['Hsp_query-to']:
            return None, None

        align_index = pos - self._hsp['Hsp_query-from']
        res_temp = self._hsp['Hsp_hseq'][align_index]

        # res_temp is gap, return None, None
        if not re.match('[A-Za-z]', res_temp):
            return None, None
        pos_temp = self._hsp['hit_seq_map'][align_index] + self._hsp['Hsp_hit-from']

        return res_temp, pos_temp


    def write(self, fname, file_format='txt', mode='w'):
        '''
        Write alignment data to file or output stream.

        Parameters
        ----------
        fname: str
            File name to write to.
        file_format: str
            One of 'txt' or 'xml'.
        mode: str
            Specify write mode. One of 'w', 'a'.
        '''

        # Check arguments
        if mode not in ['a', 'w']:
            raise RuntimeError('{} is an invalid mode.'.format(mode))

        if self._empty:
            return

        with open(fname, mode) as outF:
            if file_format == 'txt':
                # print header
                raise NotImplementedError('txt method not implemented yet...\nUse xml instead.')
            elif file_format == 'xml':
                outF.write(ET.tostring(self._tree, encoding='unicode'))
                if mode == 'a':
                    outF.write('\n')
            else:
                raise RuntimeError('{} is an unknown file_format!'.format(file_format))


def _blastp_worker(search_item, db = None, verbose=False):

    #query = '>sp|{}|{}\n{}'.format(search_item[0], search_item[2], search_item[3])
    query = search_item[3]
    return_code, dat = blastp(search_item[1], db, query, verbose=verbose)
    return dat


def align_all(peptides, sequences, db_path, organisms, nThread=None, show_bar=True, verbose=False):

    #construct list to pass to blastp worker
    search_list = list()
    for id in set([x['id'] for x in peptides]):
        for o in organisms:
            # search_list is tuple of (id, organisms, description, sequence)
            search_list.append((id, o, sequences[id][0], sequences[id][1]))

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
        length = len(search_list)
        for i, it in enumerate(search_list):
            sys.stdout.write('Working on {} of {}'.format(i, length))
            results.append(_blastp_worker(it, db=db_path, verbose=verbose))

    assert(len(search_list) == len(results))

    ret=dict()
    for sl, r in zip(search_list, results):
        if sl[0] not in ret:
            ret[sl[0]]=dict()
        ret[sl[0]][sl[1]]=Alignment(r, query_id = sl[0], query_description = sl[2], query_organism = sl[1])

    return ret

