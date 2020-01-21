
import argparse

PARENT_PARSER = argparse.ArgumentParser(add_help=False)

PARENT_PARSER.add_argument('-f', '--file_type', default = 'cimage', choices=['cimage','dtaselect'],
                    help='Choose input file format. cimage is the default.')

PARENT_PARSER.add_argument('-s', '--write_seq', choices=[0,1], type=int, default=0,
                    help='Write protein sequences in input to fasta file? 0 is the default.')

PARENT_PARSER.add_argument('--ofname', default='Cysteine_annotation.tsv',
                    help='Name of file to write results to.')

PARENT_PARSER.add_argument('-a', '--align', choices=[0,1], type=int, default=0,
                    help='Choose whether to balast protein sequences to determine cysteine conservation. '
                         'If this option is specified, a database dir must also be specified with the --database_dir option. '
                         '0 is the default.')

PARENT_PARSER.add_argument('-w', '--write_alignment_data', choices=[0,1], type=int, default=0,
                    help='Choose whether to write alignment data. '
                         '0 is the default.')

PARENT_PARSER.add_argument('-d', '--database_dir', type = str,
                    help = 'Path to directory containing sequence databases to use for alignment.')

PARENT_PARSER.add_argument('--residue_sep', type = str, default='|',
                    help = 'Residue seperator for peptides with multiple modified residues. Default is \'|\'.')

PARENT_PARSER.add_argument('--fxn_sep', type = str, default='!',
                    help = 'Function seperator for peptides with multiple modified residues. Default is \'!\'.')

PARENT_PARSER.add_argument('-o', '--defined_organism', default='none', type=str)

PARENT_PARSER.add_argument('-v', '--verbose', choices=[0,1], type=int, default=0,
                    help='Print verbose output?')

PARENT_PARSER.add_argument('input_file', type = str, help = 'Path to input file.')

