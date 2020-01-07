
import re


def write_fasta_entry(fname, acession, sequence, description = '', append = True):
    '''
    Write fasta entry to `fname`

    Paramaters
    ----------
    fname: str
        Path to file to write to.
    acession: str
        Unique identifer for sequence.
    sequence: str
        Entry sequence.
    description: str
        Entry description (optional)
    append: bool
        Should `fname` be overwitten or appended to?
    '''

    with open(fname, 'a' if append else 'w') as outF:
        outF.write('\n>sp|{}|{}\n{}'.format(acession, description, sequence))


class FastaFile(object):
    '''
    Basic FastaFile container.

    A file is parsed with a regular expression for fasta entries.
    (see self.entry_re for the regex wich is used.)
    The text of the file is stored in self._fbuff, and
    a dict containing ids, and index offsets is stored in self._id_offsets.
    '''

    uniprot_id_re = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
    sequence_re = '[A-Za-z\s]+'
    entry_re = '>[st][rp]\|({}\|.*\n{})[>]?'.format(uniprot_id_re, sequence_re)

    def __init__():
        self._id_offsets = dict()
        self._fbuff = str()

    def _find_offsets(self):
        for m in re.finditer(self.entry_re, self._fbuff):
            self._id_offsets[m.group(1)] = m.span()

    def read(self, fname):
        '''
        Read `fname` and populate self._id_offsets
        '''

        with open(fname, 'r') as inF:
            self._fbuff = inF.read()
        self._find_offsets()

    def id_exists(self, acession):
        '''
        Return True if `acession` exists file.
        '''

        return acession in self._id_offsets

    def get_offset(acession):
        '''
        Get the index offsets of `acession` in self._fbuff.

        Return
        ------
        offset: tuple
            Tuple with the format (<begin>, <end>)

        Raises
        ------
        AssertionError if !self.id_exists(acession)
        '''

        assert(self.id_exists(acession))
        return self._id_offsets(acession)

    def get_sequence(self, acession):
        '''
        Returns sequence of `acession`.

        Raises
        ------
        AssertionError if !self.id_exists(acession)
        AssertionError if not able to match self.entry_re at offset.
        AssertionError if acession of match does not match `acession`
        '''

        _offset = self.get_offset(acession)
        m = re.find(self.entry_re, self._fbuff[_offset[0]:_offset[1]])
        assert(m)
        assert(m.group(1) == acession)
        return re.sub('\s', '', m.group(3))


