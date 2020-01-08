#! /usr/bin/env python

import os
import os.path
import sys
import argparse

from .submodules import MSParser, UniProt, Blast, Alignments

# List of organisms for conservation analysis
organism_list = ['human', 'mouse', 'fly', 'yeast', 'mustard', 'worms']


def parse_input(fname, file_type, defined_organism):
    _parser = None
    if file_type == 'cimage':
        _parser = MSParser.cimage
    elif file_type == 'dtaselect':
        _parser = MSParser.dtaselect
    else:
        raise RuntimeError('{} is an unknown input file_type'.format(file_type))

    with open(fname, 'r') as inF:
        lines = inF.readlines()

    ret = list()
    for line in lines:
        ret.append(_parser(line, organism_list, defined_organism))

    return ret


def main():
    parser = argparse.ArgumentParser(prog = 'cimage_annotation',
                                     description = 'Annotate functional cysteine residues in cimage output.')

    parser.add_argument('-f', '--file_type', default = 'cimage', choices=['cimage','dtaselect'],
                        help='Choose input file format. cimage is the default.')

    parser.add_argument('-s', '--write_seq', default = False, action='store_true',
                        help='Write protein sequences in input to fasta file?')

    parser.add_argument('--align', choices=[0,1], type=int, default=0,
                        help='Choose whether to balast protein sequences to determine cysteine conservation.'
                             '0 is the default.')

    parser.add_argument('-d', '--database_dir', type = str,
                        help = 'Path to directory containing sequence databases to use for allignment.')

    parser.add_argument('-o', '--defined_organism', default='none')

    parser.add_argument('-p', '--parallel', choices=[0,1], type=int, default=1,
                        help='Choose whether internet queries and protein alignments should be performed in parallel.'
                        'Parallel processing is performed on up to the number of logical cores on your system.'
                        '1 is the default.')

    parser.add_argument('input_file', type = str, help = 'Path to input file.')

    args = parser.parse_args()

    # Open cimage or dtaselect file
    file_input = parse_input(args.input_file, args.file_type, args.defined_organism)

    # Determine path for file output
    path = '{}/'.format(os.path.abspath(os.path.dirname(args.input_file)))

    header = []
    peptides = []
    proteins = []
    index = ''
    uniuque_ids = set()
    peptide_indecies = list()
    sequences = dict()

    # Analysis for a cimage file creates lists (above) of dictionaries
    for i, line in enumerate(file_input):

        if line['index'].strip() == 'index':  # create header line(s)
            header.append(line)

        elif line['index'].strip() != '':  # create peptide lines (no protein information - includes overall ratio)
            index = line['index'].strip()
            peptides.append(line)

        else:
            file_input[i]['index'] = index # create protein lines
            proteins.append(line)
            uniuque_ids.add(line['id'])
            peptide_indecies.append(i)

    sys.stdout.write('Retreiving protein Uniprot records.')
    record_dict = UniProt.get_uniprot_records(uniuque_ids, args.parallel)

    seq_written=False
    for i in peptide_indecies:
        #sys.stdout.write('Working on {}\n'.format(file_input[i]['id']))
        # Get and Parse Uniprot entry for protein
        UniProt_data = UniProt.ExPasy(file_input[i]['id'],
                                      file_input[i]['sequence'].split('.')[1].split('*')[0] +
                                      file_input[i]['sequence'].split('.')[1].split('*')[1],
                                      file_input[i]['sequence'].split('.')[1].find('*'),
                                      record_dict[file_input[i]['id']])

        file_input[i]['position'] = UniProt_data[1]   # cysteine position
        file_input[i]['cys_function'] = UniProt_data[2] # cysteine function (if known)
        file_input[i]['protein_location'] = UniProt_data[4] # protein subcellular localization (if known)

        if args.write_seq:
            if file_input[i]['id'] not in sequences: #only write sequence if it is not currently in file.
                Fasta.write_fasta_entry(SEQ_PATH,
                                        file_input[i]['id'],
                                        UniProt_data[3],
                                        description=file_input[i]['description'],
                                        append=seq_written)
                seq_written = True
        sequences[file_input[i]['id']] = UniProt_data[3]

    if args.align:
        for i in peptide_indecies:

            sys.stdout.write('Working on {}\n'.format(file_input[i]['id']))

            # blast protein sequence against each fasta database and parse those alignments to determine cysteine conservation
            for organism in organism_list:
                Blast.blastp(organism, path, args.database_dir, sequences[file_input[i]['id']])
                Alignments.align_write(path, organism)
                alignment_values = Alignments.Align_file_parse(file_input[i]['position'], path)
                file_input[i][str(organism) + '_conserved'] = alignment_values[2]

                # for comparative organism analyze Uniprot entry of best blast hit
                if organism == args.defined_organism.lower():
                    if alignment_values[0] != '':
                        if alignment_values[3] == '':
                            position = 0
                        else:
                            position = alignment_values[3]
                        Organism_data = UniProt.ExPasy_alt(alignment_values[0], position)
                        file_input[i][args.defined_organism + '_id'] = alignment_values[0]
                        file_input[i][args.defined_organism + '_evalue'] = alignment_values[1]
                        file_input[i][args.defined_organism + '_description'] = Organism_data[0]
                        file_input[i][args.defined_organism + '_position'] = str(alignment_values[3])
                        file_input[i][args.defined_organism + '_function'] = Organism_data[1]
                    else:
                        file_input[i][args.defined_organism + '_id'] = ''
                        file_input[i][args.defined_organism + '_evalue'] = ''
                        file_input[i][args.defined_organism + '_description'] = ''
                        file_input[i][args.defined_organism + '_position'] = ''
                        file_input[i][args.defined_organism + '_function'] = ''
                else:
                    pass

    # file output
    output = open(path + 'Cysteine_annotation.tsv', 'w')
    output.write(MSParser.output(header[0], organism_list, args.defined_organism))
    for peptide in peptides:
        output.write(MSParser.output(peptide, organism_list, args.defined_organism))
        for protein in proteins:
            if protein['index'].strip() == peptide['index'].strip():
                output.write(MSParser.output(protein, organism_list, args.defined_organism))
            else:
                pass
    output.close()

if __name__ == '__main__':
    main()

