
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

    acession_re = r'\w+'
    sequence_re = r'[A-Za-z\s]*'
    entry_re = r'>[st][rp]\|({})\|.*\n({})(?=[>])?'.format(acession_re, sequence_re)

    def __init__(self):
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

    def iter_ids(self):
        l = len(self._id_offsets)
        for i, k in enumerate(self._id_offsets.keys()):
            if i < l:
                yield k
            else:
                raise StopIteration

    def iter_items(self):
        l = len(self._id_offsets)
        for i, v in enumerate(self._id_offsets.items()):
            if i < l:
                yield v
            else:
                raise StopIteration

    def id_exists(self, acession):
        '''
        Return True if `acession` exists file.
        '''

        return acession in self._id_offsets

    def get_offset(self, acession):
        '''
        Get the index offsets of `acession` in self._fbuff.

        Return
        ------
        offset: tuple
            Tuple with the format (<begin>, <end>)

        Raises
        ------
        KeyError if !self.id_exists(acession)
        '''

        if not self.id_exists(acession):
            raise KeyError('{} does not exist in FastaFile!'.format(acession))
        return self._id_offsets[acession]

    def get_sequence(self, acession):
        '''
        Returns sequence of `acession`.

        Raises
        ------
        KeyError if !self.id_exists(acession)
        AssertionError if not able to match self.entry_re at offset.
        AssertionError if acession of match does not match `acession`
        '''

        _offset = self.get_offset(acession)
        m = re.search(self.entry_re, self._fbuff[_offset[0]:_offset[1]])
        assert(m)
        assert(m.group(1) == acession)
        return re.sub('\s', '', m.group(2))


