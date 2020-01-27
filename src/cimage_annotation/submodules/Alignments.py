
import sys
import re
import xml.etree.ElementTree as ET
from multiprocessing import Pool
from multiprocessing import cpu_count
import itertools
from tqdm import tqdm
import functools
from math import ceil

from .Blast import blastp


class Alignment(object):
    '''
    Read BLAST XML file and provide methods to access underlying alignment data.
    '''

    _XML_HEADER_ELEMENTS = {'query_length': './BlastOutput_iterations/Iteration/Iteration_query-len',
                           'query_description': './BlastOutput_iterations/Iteration/Iteration_query-def',
                           'query_id': './BlastOutput_iterations/Iteration/Iteration_query-ID'}

    _XML_HITS_PATH = './BlastOutput_iterations/Iteration/Iteration_hits/Hit'
    _XML_HSP_PATH = 'Hit_hsps/Hsp/'
    _XML_QUERY_SEQ_NAME = 'Hsp_qseq'
    _XML_HIT_SEQ_NAME = 'Hsp_hseq'
    _XML_MIDLINE_SEQ_NAME = 'Hsp_midline'
    _XML_QUERY_SEQ_PATH = _XML_HSP_PATH + 'Hsp_qseq'
    _XML_HIT_SEQ_PATH = _XML_HSP_PATH + 'Hsp_hseq'
    _XML_MIDLINE_SEQ_PATH = _XML_HSP_PATH + 'Hsp_midline'
    _ALIGNMENT_LINE_LENGTH = 60
    _VERBOSE = False

    _XML_MATCH_ELEMENTS = {'match_num': 'Hit_num',
                          'query_from': 'Hit_hsps/Hsp/Hsp_query-from',
                          'query_to': 'Hit_hsps/Hsp/Hsp_query-to',
                          'match_id': 'Hit_accession',
                          'match_description': 'Hit_def',
                          'match_evalue': 'Hit_hsps/Hsp/Hsp_evalue',
                          'match_from': 'Hit_hsps/Hsp/Hsp_hit-from',
                          'match_to': 'Hit_hsps/Hsp/Hsp_hit-to',
                          'match_length': 'Hit_hsps/Hsp/Hsp_align-len'}


    def __init__(self, raw_xml, query_id=None, query_description=None, query_organism=None) :
        '''
        Default constructor

        Parameters
        ----------
        raw_xml: str
            BLAST XML file as text.
        query_id: str
            Acession of query.
        query_description: str
            Description of query.
        query_organism: str
            Query organism.
        '''

        self.query_id = query_id
        self.query_description = query_description
        self.query_organism = query_organism

        if raw_xml:
            self._tree = ET.fromstring(raw_xml)
            self._best_hit = self._tree.find(self._XML_HITS_PATH)
            self._hsp = None
            self._add_text_to_path('BlastOutput_iterations/Iteration/Iteration_query-ID',
                                   self.query_id,
                                   error_level=1 if self._VERBOSE else 0)
            self._add_text_to_path('BlastOutput_iterations/Iteration/Iteration_query-def',
                                   self.query_description,
                                   error_level=1 if self._VERBOSE else 0)

            if self._best_hit is None:
                self._empty = True
            else:
                self._empty = False
        else:
            self._empty = True


    def _add_text_to_path(self, path, text, error_level = 0):
        '''
        Add path to element in self._tree

        Parameters
        ----------
        path: str
            Xpath to element
        text: str
            Test to add to element
        error_level: int
            Specify what to do if `path` does not exist.
            0: do nothing, exit from function silently.
            1: Print a warning and exit.
            2: Raise an AttributeError.
        '''

        try:
            self._tree.find(path).text = text
        except AttributeError as e:
            if error_level != 0:
                sys.stderr.write('WARN: No element at path {}'.format(path))
            if error_level != 1 or error_level != 0:
                raise AttributeError(e)


    def _get_hit_path_text(self, path):
        temp = self._best_hit.find(path)
        return None if temp is None else temp.text


    def get_best_id(self):
        '''
        Get acession of best alignment match.

        Returns
        -------
        best_id: str
        '''
        if self._empty:
            return ''
        return self._get_hit_path_text('./Hit_accession')


    def get_best_description(self):
        '''
        Get description of best alignment match.

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
        self._hsp = {x.tag: x.text for x in self._best_hit.findall(self._XML_HSP_PATH)}

        self._hsp['hit_seq_map'] = list()
        seq_index = 0
        for i, c in enumerate(self._hsp[self._XML_HIT_SEQ_NAME]):
            self._hsp['hit_seq_map'].append(seq_index)
            if re.match('[A-Za-z]', c):
                seq_index += 1

        self._hsp['query_seq_map'] = list()
        for i, c in enumerate(self._hsp[self._XML_QUERY_SEQ_NAME]):
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
        query_residue = self._hsp[self._XML_QUERY_SEQ_NAME][query_index]
        hit_residue = self._hsp[self._XML_HIT_SEQ_NAME][query_index]
        if self._VERBOSE:
            sys.stdout.write('{}: {} == {} -> {}\n'.format(pos,
                                                           query_residue,
                                                           hit_residue,
                                                           (hit_residue == query_residue)))

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
        res_temp = self._hsp[self._XML_HIT_SEQ_NAME][align_index]

        # res_temp is gap, return None, None
        if not re.match('[A-Za-z]', res_temp):
            return None, None
        pos_temp = self._hsp['hit_seq_map'][align_index] + self._hsp['Hsp_hit-from']

        return res_temp, pos_temp


    @staticmethod
    def _write_element(out, tag, name, value):
        out.write('{}\t{}\t{}\n'.format(tag, name, value))


    @staticmethod
    def _search_path_text(result, empty_value = ''):
        return empty_value if result is None else result.text


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
            if self._VERBOSE:
                sys.stderr.write('WARN: Attempting to write an empty Alignment file.')
            return

        with open(fname, mode) as outF:
            if file_format == 'txt':
                # print header
                for name, path in self._XML_HEADER_ELEMENTS.items():
                    text = self._search_path_text(self._tree.find(path))
                    if text == '' and self._VERBOSE:
                        sys.stderr.write('WARN: No element at path {}'.format(path))
                    self._write_element(outF, 'H', name, text)

                # print hits
                for hit in self._tree.findall(self._XML_HITS_PATH):
                    outF.write('\n')
                    # print hit header
                    for name, path in self._XML_MATCH_ELEMENTS.items():
                        text = self._search_path_text(hit.find(path))
                        if text == '' and self._VERBOSE:
                            sys.stderr.write('WARN: No element at path {}'.format(path))
                        self._write_element(outF, 'M', name, text)

                    # print alignment data
                    hit_seq = self._search_path_text(hit.find(self._XML_HIT_SEQ_PATH))
                    query_seq = self._search_path_text(hit.find(self._XML_QUERY_SEQ_PATH))
                    midline_seq = self._search_path_text(hit.find(self._XML_MIDLINE_SEQ_PATH))
                    max_name_len = max([len(x) for x in [self._XML_HIT_SEQ_NAME, self._XML_QUERY_SEQ_NAME, self._XML_MIDLINE_SEQ_NAME]])
                    length = len(query_seq)
                    n_lines = ceil(length / self._ALIGNMENT_LINE_LENGTH)
                    for i in range(n_lines):
                        begin = i * self._ALIGNMENT_LINE_LENGTH
                        end = (i + 1) * self._ALIGNMENT_LINE_LENGTH
                        end = end if end < length else length

                        self._write_element(outF, 'A', '{}{}'.format(self._XML_QUERY_SEQ_NAME,
                                                                     ' ' * (max_name_len - len(self._XML_QUERY_SEQ_NAME))),
                                            '{} {}'.format(query_seq[begin:end], end))

                        self._write_element(outF, 'A', '{}{}'.format(self._XML_MIDLINE_SEQ_NAME,
                                                                     ' ' * (max_name_len - len(self._XML_MIDLINE_SEQ_NAME))),
                                            '{} {}'.format(midline_seq[begin:end], end))

                        self._write_element(outF, 'A', '{}{}'.format(self._XML_HIT_SEQ_NAME,
                                                                     ' ' * (max_name_len - len(self._XML_HIT_SEQ_NAME))),
                                            '{} {}'.format(hit_seq[begin:end], end))
                        outF.write('\n')

            elif file_format == 'xml':
                outF.write(ET.totexting(self._tree, encoding='unicode'))
                if mode == 'a':
                    outF.write('\n')
            else:
                raise RuntimeError('{} is an unknown file_format!'.format(file_format))


def _blastp_worker(search_item, db = None, verbose=False):
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
    Alignment._VERBOSE = verbose
    for sl, r in zip(search_list, results):
        if sl[0] not in ret:
            ret[sl[0]]=dict()
        ret[sl[0]][sl[1]]=Alignment(r, query_id=sl[0],
                                    query_description=sl[2],
                                    query_organism=sl[1])

    return ret

