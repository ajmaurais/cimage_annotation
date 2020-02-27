
import sys
import csv

class DataFrame(object):
    '''
    Lightweight DataFrame class which recreates some of the
    functionality as a Pandas.DataFrame without having to
    import the entire Pandas library.
    '''

    def __init__(self, data=None):
        self.nrow = 0
        self.ncol = 0
        self._data = list()
        self.columns = list()
        self._keys = dict()

    def empty(self):
        ''' Return true if DataFrame is empty. '''
        return len(self._data) == 0

    def __len__(self):
        return self.nrow

    def __getitem__(self, col):
        if col not in self._keys:
            raise KeyError('{} is not a column in DataFrame!'.format(col))
        return self._data[self._keys[col]]

    def __setitem__(self, col, value):
        if not isinstance(value, list):
            raise ValueError('value must be if type {}. Recieved {}'.format(type(self._data), type(value)))
        if not self.empty():
            if len(value) != len(self._data[0]):
                raise ValueError('Atempting to add column of lenght {}'
                    ' to DataFrame with columns of length {}'.format(len(value), len(self._data[0])))
        if col not in self._keys:
            self.columns.append(col)
            self._keys[col] = len(self.columns) - 1
            self.ncol += 1
            self._data.append(list())
        self._data[self._keys[col]] = value

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
                yield i, {col:self.__getitem__(col)[i] for col in self.columns}

    def to_csv(self, fname, sep='\t'):
        '''
        Write DataFrame to file.

        Paramaters
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
    File format and deliminators are detected automatically using csv.Sniffer().

    Parameters
    ----------
    fname: str
        Path to file to read.
    hasHeader: bool
        Does the first line in the file contain column headers?
    delim: str
        Column deliminator.

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

    #add row keys
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
        _keys = [x for x in range(len(lines[0]))]
        _start = 0
    ret._data = [list() for x in _keys]
    ret.columns = _keys
    ret._keys = {k: i for i, k in enumerate(_keys)}
    ret.ncol = len(_keys)

    #iterate through lines
    for i, line in enumerate(lines[_start:]):
        for j, elem in enumerate(line):
            ret._data[j].append(elem)
    ret.nrow = len(ret._data[0])

    return ret

