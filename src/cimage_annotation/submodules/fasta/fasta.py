
import re


def write_fasta_entry(fname, acession, sequence, description='', write_mode='a'):
    '''
    Write fasta entry to `fname`

    Parameters
    ----------
    fname: str
        Path to file to write to.
    accession: str
        Unique identifier for sequence.
    sequence: str
        Entry sequence.
    description: str
        Entry description (optional)
    write_mode: str
        File write mode (must be 'w' or 'a')

    Raises
    ------
    FileNotFoundError:
        If file can not be written to.
    ValueError:
        If invalid write mode.
    '''

    if write_mode not in ('a', 'w'):
        raise ValueError('{} is an invalid write_mode.'.format(write_mode))
    with open(fname, write_mode) as outF:
        outF.write('\n>sp|{}|{}\n{}'.format(acession, description, sequence))


class FastaFile(object):
    '''
    Basic FastaFile container, optimized for random access.

    A file is parsed with a regular expression for UniProt fasta entries.
    (see FastaFile.entry_re for the regex which is used.)
    The text of the file is stored in FastaFile._fbuff, and a dict
    containing ids, and index offsets is stored in FastaFile._id_offsets.

    Examples
    --------
    Initialize and read FastaFile.

    >>> fasta = FastaFile()
    >>> fasta.read(file_path)

    You can now acess the sequences in the fasta file

    >>> fasta.get_sequence('P26641')
    'MAAGTLYTYPENWRAFKALIAAQYSGAQVRVLSAPPHFHFGQTNRTPEFLRKFPAGKVPAFEGDDGFCVFESNAI
    AYYVSNEELRGSTPEAAAQVVQWVSFADSDIVPPASTWVFPTLGIMHHNKQATENAKEEVRRILGLLDAYLKTRTF
    LVGERVTLADITVVCTLLWLYKQVLEPSFRQAFPNTNRWFLTCINQPQFRAVLGEVKLCEKMAQFDAKKFAETQPK
    KDTPRKEKGSREEKQKPQAERKEEKKAAAPAPEEEMDECEQALAAEPKAKDPFAHLPKSTFVLDEFKRKYSNEDTL
    SVALPYFWEHFDKDGWSLWYSEYRFPEELTQTFMSCNLITGMFQRLDKLRKNAFASVILFGTNNSSSISGVWVFRG
    QELAFPLSPDWQVDYESYTWRKLDPGSEETQTLVREYFSWEGAFQHVGKAFNQGKIFK'

    Attempting to acess a sequence which doesn't exist will raise a KeyError

    >>> fasta.get_sequence('DUMMY_ID')
    KeyError: 'DUMMY_ID does not exist in FastaFile!'

    You can check wether a sequence exists in the file for a given acession.

    >>> fasta.id_exists('DUMMY_ID')
    False
    >>> fasta.id_exists('P26641')
    True

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
        '''
        Iterate over FastaFile entry IDs.

        Yields
        ------
        index: int
            Index of element
        id: str
            Accession of entry.
        '''

        for i, k in enumerate(self._id_offsets.keys()):
            yield i, k

    def iter_items(self):
        '''
        Iterate over FastaFile entries.

        Yields
        ------
        index: int
            Index of element
        id: str
            Accession of entry.
        seq: str
            Sequence of entry.
        '''

        for i, (k, v) in enumerate(self._id_offsets.items()):
            yield i, k, self._get_sequence(k, v[0], v[1])

    def __contains__(self, acession):
        return acession in self._id_offsets

    def id_exists(self, acession):
        '''
        Return True if `accession` exists file.
        '''
        return self.__contains__(acession)

    def _get_offset(self, acession):
        '''
        Get the index offsets of `accession` in self._fbuff.

        Return
        ------
        offset: tuple
            Tuple with the format (<begin>, <end>)

        Raises
        ------
        KeyError if !self.id_exists(accession)
        '''

        if not self.id_exists(acession):
            raise KeyError('{} does not exist in FastaFile!'.format(acession))
        return self._id_offsets[acession]
    
    def _get_sequence(self, acession, begin, end):
        m = re.search(self.entry_re, self._fbuff[begin:end])
        assert(m)
        assert(m.group(1) == acession)
        return re.sub('\s', '', m.group(2))

    def get_sequence(self, acession):
        '''
        Returns sequence of `accession`.

        Raises
        ------
        KeyError:
            if !self.id_exists(accession)
        AssertionError:
            if not able to match self.entry_re at offset.
        AssertionError:
            if accession of match does not match `accession`
        '''

        _offset = self._get_offset(acession)
        return self._get_sequence(acession, _offset[0], _offset[1])

