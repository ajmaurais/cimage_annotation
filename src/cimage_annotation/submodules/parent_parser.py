
import argparse

PARENT_PARSER = argparse.ArgumentParser(add_help=False)

PARENT_PARSER.add_argument('-f', '--file_type', default='cimage', choices=['cimage', 'tsv'],
                           help='Choose input file format. cimage is the default.'
                                 'For tsv format, there must be columns for protein ID and peptide sequence.')

PARENT_PARSER.add_argument('--id_col', default='id',
                           help='Specify the columm header containing protein Uniprot IDs for tsv input.')

PARENT_PARSER.add_argument('--seq_col', default='sequence',
                           help='Specify the columm header containing peptide sequences for tsv input.')

PARENT_PARSER.add_argument('--description_col', default='description',
                           help='Specify the columm header containing protein descriptions for tsv input. '
                                'Only required for query description in alignment output files.')

PARENT_PARSER.add_argument('-s', '--write_seq', action='store_true', default=False,
                           help='Write protein sequences in input to fasta file? 0 is the default.')

PARENT_PARSER.add_argument('--all_features', action='store_true', default=False,
                           help='Should all residue feature annotations be included? '
                                'By default, only a simplified set of features are included.')

PARENT_PARSER.add_argument('--ofname', default='residue_annotation.tsv',
                           help='Name of file to write results to.')

PARENT_PARSER.add_argument('-a', '--align', action='store_true', default=False,
                           help='Choose whether to balast protein sequences to determine residue conservation. '
                                'If this option is specified, a database dir must also be specified with the --database_dir option. '
                                '0 is the default.')

PARENT_PARSER.add_argument('-w', '--write_alignment_data', action='store_true', default=False,
                           help='Choose whether to write alignment data. '
                                '0 is the default.')

PARENT_PARSER.add_argument('--align_format', choices=['xml', 'txt'], default='txt',
                           help='Alignment file format. xml is BLAST XML format, sutible for programming, '
                           'txt is human readable.')

PARENT_PARSER.add_argument('--evalue_co', default=1e-5, type=float,
                           help='Alignment evalue cuttoff for a residue to be considered conserved. 1e-5 is the default.')

PARENT_PARSER.add_argument('-d', '--database_dir', type=str,
                           help='Path to directory containing sequence databases to use for alignment.')

PARENT_PARSER.add_argument('-o', '--defined_organism', default='none', type=str,
                           help='Define organism to look up function of conserved residues.')

PARENT_PARSER.add_argument('-v', '--verbose', action='store_true', default=False,
                           help='Print verbose output?')

PARENT_PARSER.add_argument('input_file', type=str, help='Path to input file.')

