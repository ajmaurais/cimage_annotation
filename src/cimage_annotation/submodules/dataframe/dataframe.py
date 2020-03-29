
import sys
import csv


class DataFrame(object):
    '''
    Lightweight DataFrame class which recreates some of the
    functionality as a Pandas.DataFrame without having to
    import the entire Pandas library.

    Parameters
    ----------
    data: dict, default None
        Dict with column names as keys, and column values as list like objects.

    Examples
    --------
    Read .tsv file into DataFrame

    >>> dat = dataframe.read_tsv('testData/input_peptides.tsv')
    >>> dat
                    ID                      seq
    0          'P26641' 'FPEELTQTFMSC*NLITGMFQR'
    1          'Q9NTZ6'        'VC*AHITNIPFSITK'
    2          'P19012'         'AGLENSLAETEC*R'
    3          'P42224'     'HLLPLWNDGC*IMGFISK'
    4          'Q15257'        'QSVSC*DEC*IPLPR'
    ########## Reached self._MAX_REPR_PRINT ##########
    
    Get a list of column names by calling the columns member.
    
    >>> dat.columns
    ['ID', 'seq']

    Get and set column contents with the [] operator.

    >>> dat['ID']
    ['P26641',
     'Q9NTZ6',
     'P19012',
     'P42224',
     'Q15257',
        ...
    >>> dat['another_column'] = ['' for _ in range(dat.nrow)]
    >>> dat
           ID                      seq another_column
    0 'P26641' 'FPEELTQTFMSC*NLITGMFQR'             ''
    1 'Q9NTZ6'        'VC*AHITNIPFSITK'             ''
    2 'P19012'         'AGLENSLAETEC*R'             ''
    3 'P42224'     'HLLPLWNDGC*IMGFISK'             ''
    4 'Q15257'        'QSVSC*DEC*IPLPR'             ''
    ########## Reached self._MAX_REPR_PRINT ##########

    Iterate through the DataFrame contents by row with the iterrows method.

    >>> for i, r in dat.iterrows():
    ...    print(r['seq'])
    FPEELTQTFMSC*NLITGMFQR
    VC*AHITNIPFSITK
    AGLENSLAETEC*R
    HLLPLWNDGC*IMGFISK
    QSVSC*DEC*IPLPR
    DFTPVC*TTELGR
        ...
    '''

    _MAX_REPR_PRINT = 100

    def __init__(self, data=None):
        self.columns = list()
        self._keys = dict()
        self._data = list()
        self.nrow = 0
        self.ncol = 0
        if data is not None:
            lengths = set()
            for i, (k, v) in enumerate(data.items()):
                lengths.add(len(v))
                if len(lengths) != 1:
                    raise ValueError('arrays must all be same length')
                self.columns.append(k)
                self._keys[k] = i
                self._data.append(v)
            self.nrow = len(self._data[0])
            self.ncol = len(self.columns)

    def empty(self):
        ''' Return true if DataFrame is empty. '''
        return len(self._data) == 0

    def __len__(self):
        return self.nrow

    def __getitem__(self, col):
        if isinstance(col, list):
            ret = DataFrame()
            for c in col:
                if c not in self._keys:
                    raise KeyError('{} is not a column in DataFrame!'.format(col))
                ret[c] = self._data[self._keys[c]]
        else:
            if col not in self._keys:
                raise KeyError('{} is not a column in DataFrame!'.format(col))
            ret = self._data[self._keys[col]]
        return ret

    def __setitem__(self, col, value):
        if not isinstance(value, list):
            raise ValueError('value must be of type {}. Received {}'.format(type(self._data),
                                                                            type(value)))
        if not self.empty():
            if len(value) != len(self._data[0]):
                raise ValueError('Attempting to add column of length {}'
                                 ' to DataFrame with columns of length {}'.format(len(value),
                                                                                  len(self._data[0])))
        else:
            self.nrow = len(value)

        if col not in self._keys:
            self.columns.append(col)
            self._keys[col] = len(self.columns) - 1
            self.ncol += 1
            self._data.append(list())
        self._data[self._keys[col]] = value

    def _format_r_rows(self, r, print_header=True, delim=' '):
        '''
        Return formatted string with rows in range `r`.
        '''
        
        spaces = list()
        for c in self.columns:
            col_index = self._keys[c]
            max_temp = max([len(repr(self._data[col_index][i])) for i in r])
            spaces.append(max(len(c), max_temp))

        index_space = max(len(str(x)) for x in r)

        ret = ''
        if print_header:
            ret += '{}{}\n'.format(' ' * index_space,
                    delim.join([c.rjust(spaces[c_i]) for c_i, c in enumerate(self.columns)]))
        for i in r:
            ret += '{}{}{}{}\n'.format(i, ' ' * (index_space - len(str(i))), delim,
                    delim.join([repr(self._data[self._keys[x]][i]).rjust(spaces[c_i]) for c_i, x in enumerate(self.columns)]))
        return ret

    def __repr__(self):
        print_n = self._MAX_REPR_PRINT if self._MAX_REPR_PRINT < self.nrow else self.nrow
        ret = ''
        ret += self._format_r_rows(range(print_n), sys.stdout)
        if print_n < self.nrow:
            ret += '{0} Reached self._MAX_REPR_PRINT {0}'.format('#' * 10)
        return ret

    def head(self, n = 5, delim=' '):
        '''
        Print the first `n` rows.
        '''
        sys.stdout.write(self._format_r_rows(range(n), delim=delim))


    def iterrows(self):
        '''
        Iterate over DataFrame rows as (index, dict) pairs.

        Yields
        ------
        index: int
            Row number starting from 0.
        dat: dict
            Dictionary with data from row.
        '''
        for i in range(self.nrow):
            if i < self.nrow:
                yield i, {col: self.__getitem__(col)[i] for col in self.columns}

    def to_tsv(self, fname, sep='\t'):
        '''
        Write DataFrame to file.

        Parameters
        ----------
        fname: str
            Path to file to write to.

        Raises
        ------
        FileNotFoundError
            If output directory does not exist.
        '''

        with open(fname, 'w') as outF:
            outF.write('{}'.format(sep).join(self.columns))
            outF.write('\n')
            for i in range(self.nrow):
                for row_i, c in enumerate(self.columns):
                    if row_i == 0:
                        outF.write(self.__getitem__(c)[i])
                    else:
                        outF.write('{}{}'.format(sep, self.__getitem__(c)[i]))
                outF.write('\n')


def read_tsv(fname, hasHeader=True):
    '''
    Read csv or tsv file into DataFrame.
    File format and delimiters are detected automatically using csv.Sniffer().

    Parameters
    ----------
    fname: str
        Path to file to read.
    hasHeader: bool
        Does the first line in the file contain column headers?

    Return
    ------
    dataframe: DataFrame
        Object with data from `fname`.
    '''

    lines = list()
    with open(fname, 'r') as inF:
        dialect = csv.Sniffer().sniff(inF.readline())
        inF.seek(0)
        for line in csv.reader(inF, dialect):
            lines.append(line)

    ret = DataFrame()

    # add row keys
    _keys = dict()
    if hasHeader:
        _keys = [x.strip() for x in lines[0]]
        keys_temp = dict()
        for key in _keys:
            if key not in keys_temp:
                keys_temp[key] = 0
            keys_temp[key] += 1
        bad_cols = 0
        for k, v in keys_temp.items():
            if v > 1:
                sys.stderr.write('Duplicate column name found: {}\n'.format(k))
                bad_cols += 1
        if bad_cols > 0:
            raise RuntimeError('{} invalid column names in {}'.format(bad_cols, fname))
        _start = 1
    else:
        _keys = list(range(len(lines[0])))
        _start = 0
    ret._data = [list() for x in _keys]
    ret.columns = _keys
    ret._keys = {k: i for i, k in enumerate(_keys)}
    ret.ncol = len(_keys)

    # iterate through lines
    for i, line in enumerate(lines[_start:]):
        for j, elem in enumerate(line):
            ret._data[j].append(elem)
    ret.nrow = len(ret._data[0])

    return ret

