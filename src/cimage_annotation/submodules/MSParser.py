
from .Alignments import organism_list
from .dataframe import DataFrame, read_tsv

PRINT_COLS=['index', 'id', 'symbol', 'description', 'protein_location', 'sequence', 'mass', 'position', 'res_function', 'domains']
ALLIGNMENT_COLUMNS = ['id', 'evalue', 'description', 'position', 'function']
ADDED_COLUMNS = ['position', 'res_function', 'domains', 'protein_location']

'''
This file implements containers for various input file formats.
In order to be a valid input file container, the class must implement the
folowing members.

Input/output
------------
read(fname, defined_organism):
    Read input file and populate defined organism slots.
write(fname):
    write formated output file.

Data acess
----------
__len__():
    Return number of peptides (or rows).
iterpeptides():
    Iterate over peptides as (index, dict) pairs.
set_peptide_value(index, key, value):
    Set value of peptide value at index.
add_column(name):
    Add empty column with specified name.

Required members
----------------
unique_ids:
    Set containing unique protein IDs.
'''


class Cimage_file():
    '''
    Container for cimage output file.
    '''

    def __init__(self):
        self.fname = str()
        self.peptides = list()
        self.residues = list()
        self.unique_ids = set()
        self.header = dict()
        self.defined_organism = str()

    def __len__(self):
        return len(self.peptides)

    def iterpeptides(self):
        for i, p in enumerate(self.peptides):
            yield i, p

    def set_peptide_value(self, index, key, value):
        self.peptides[index][key] = value

    def add_column(self, name):
        pass

    @staticmethod
    def _parser(line, organism_list, defined_organism):
        line_list = line.split('\t')
        line_dict = {}
        dict_terms = ['index', 'id', 'description', 'symbol', 'sequence', 'mass']

        i = 0
        for item in dict_terms:
            line_dict[item] = line_list[i]
            i += 1

        n = 1
        m = len(line_dict)
        while n + m <= len(line_list):
            line_dict[n] = line_list[n + m -1]
            n += 1

        if line_dict['description'].strip() == 'description':
            line_dict['id'] = 'id'
            line_dict['protein_location'] = 'protein location'
            line_dict['position'] = 'residue position'
            line_dict['res_function'] = 'residue function'
            line_dict['domains'] = 'domains'
            line_dict[defined_organism + '_id'] = defined_organism + '_id'
            line_dict[defined_organism + '_evalue'] = defined_organism + '_evalue'
            line_dict[defined_organism + '_description'] = defined_organism + '_description'
            line_dict[defined_organism + '_position'] = defined_organism + '_position'
            line_dict[defined_organism + '_function'] = defined_organism + '_function'
            for organism in organism_list:
                line_dict[organism + '_conserved'] = organism + '_conserved'

        else:
            line_dict['protein_location'] = ''
            line_dict['position'] = ''
            line_dict['res_function'] = ''
            line_dict['domains'] = ''
            line_dict[defined_organism + '_id'] = ''
            line_dict[defined_organism + '_evalue'] = ''
            line_dict[defined_organism + '_description'] = ''
            line_dict[defined_organism + '_position'] = ''
            line_dict[defined_organism + '_function'] = ''
            for organism in organism_list:
                line_dict[organism + '_conserved'] = ''

        return line_dict

    def read(self, fname, defined_organism):
        '''
        Read and parse cimage output file.

        Paramaters
        ----------
        fname: str
            Path to file to read.
        defined_organism: str
            Defined organism.
        '''

        self.defined_organism = defined_organism
        with open(fname, 'r') as inF:
            lines = inF.readlines()

        index_temp = ''
        for i, line in enumerate(lines):
            line_dict = self._parser(line, organism_list, defined_organism)

            if line_dict['index'].strip() == 'index':  # create header line(s)
                self.header = line_dict

            elif line_dict['index'].strip() != '':  # create peptide lines (no protein information - includes overall ratio)
                index_temp = line_dict['index'].strip()
                self.residues.append(line_dict)
            else:
                self.peptides.append(line_dict)
                self.peptides[-1]['index'] = index_temp
                self.unique_ids.add(line_dict['id'])


    @staticmethod
    def _write_line(outF, line_dict, defined_organism):

        outF.write('\t'.join([line_dict[x] for x in PRINT_COLS]))

        for organism in organism_list:
            outF.write('\t{}'.format(line_dict[organism + '_conserved']))
        outF.write('\t')
        outF.write('\t'.join([str(line_dict['{}_{}'.format(defined_organism, x)]) for x in ALLIGNMENT_COLUMNS]))

        n = 1
        while n < (len(line_dict) - 20):
            outF.write(('\t{}'.format(line_dict[n].strip())))
            n += 1
        outF.write('\n')


    def write(self, fname):
        '''
        Write annotated cimage data to file.

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
            self._write_line(outF, self.header, self.defined_organism)
            for residue in self.residues:
                self._write_line(outF, residue, self.defined_organism)
                for peptide in self.peptides:
                    if peptide['index'].strip() == residue['index'].strip():
                        self._write_line(outF, peptide, self.defined_organism)

class Tsv_file():
    def __init__(self, id_col='protein_ID', seq_col='sequence'):
        self.dat = DataFrame()
        self.unique_ids = set()

        # column names
        self.id_col = id_col
        self.seq_col = seq_col

    def read(self, fname, defined_organism):
        self.dat = read_tsv(fname)
        for col in (self.id_col, self.seq_col):
            if col not in self.dat.columns:
                raise KeyError('Required column: "{}" not found!'.format(col))

        for col in ADDED_COLUMNS:
            self.add_column(col)

        self.unique_ids = set(self.dat[self.id_col])

    def __len__(self):
        return len(self.dat)

    def iterpeptides(self):
        return self.dat.iterrows()

    def write(self, fname):
        self.dat.to_tsv(fname)

    def set_peptide_value(self, index, key, value):
        self.dat[key][index] = value

    def add_column(self, name):
        self.dat[name] = ['' for _ in range(self.dat.nrow)]

class Dtaselect():
    def __init__(self):
        raise NotImplementedError()


