#! /usr/bin/env python
#
#

# Import modules
import os
import os.path
import sys
from subprograms import __MSParser__
from subprograms import __UniProt__
from subprograms import __Blast__
from subprograms import __Alignments__
# Check for all necessary sys.argv commands

# List of organisms for conservation analysis
organism_list = ['human', 'mouse', 'fly', 'yeast', 'mustard', 'worms']

# Open cimage or dtaselect file
f_input = open(sys.argv[1], 'r')
file_input = f_input.readlines()
f_input.close()

# Determine path for file output
f_path = sys.argv[1].split('/')
f_path.pop(len(f_path)-1)
path = ''
for dir in f_path:
    path += (dir + '/')

header = []
peptides = []
proteins = []
index = ''

# Analysis for a cimage file creates lists (above) of dictionaries
if sys.argv[2].lower() == 'cimage':
    for item in file_input:
        line = __MSParser__.cimage(item, organism_list, sys.argv[3])   # parse the lines of a cimage file, creating output framework (columns)

        if line['index'].strip() == 'index':  # create header line(s)
            header.append(line)

        elif line['index'].strip() != '':  # create peptide lines (no protein information - includes overall ratio)
            index = line['index'].strip()
            peptides.append(line)

        else:
            line['index'] = index # create protein lines
            proteins.append(line)
            
            sys.stdout.write('Working on {}\n'.format(line['id']))
            # Get and Parse Uniprot entry for protein
            UniProt_data = __UniProt__.ExPasy(line['id'],
                                              line['sequence'].split('.')[1].split('*')[0] +
                                              line['sequence'].split('.')[1].split('*')[1],
                                              line['sequence'].split('.')[1].find('*'))

            line['position'] = UniProt_data[1]   # cysteine position
            line['cys_function'] = UniProt_data[2] # cysteine function (if known)
            line['protein_location'] = UniProt_data[4] # protein subcellular localization (if known)

            # write full protein sequence to file for blast analysis
            f_seq = open(path + 'sequence.txt', 'w')
            f_seq.write('>sp|' + line['id'] + '|' + line['description'] + '\n' + UniProt_data[3])
            f_seq.close()

            # blast protein sequence against each fasta database and parse those alignments to determine cysteine conservation
            for organism in organism_list:
                __Blast__.blastp(organism, path)
                __Alignments__.align_write(path, organism)
                alignment_values = __Alignments__.Align_file_parse(line['position'], path)
                line[str(organism) + '_conserved'] = alignment_values[2]

                # for comparative organism analyze Uniprot entry of best blast hit
                if organism == sys.argv[3].lower():
                    if alignment_values[0] != '':
                        if alignment_values[3] == '':
                            position = 0
                        else:
                            position = alignment_values[3]
                        Organism_data = __UniProt__.ExPasy_alt(alignment_values[0], position)
                        line[sys.argv[3] + '_id'] = alignment_values[0]
                        line[sys.argv[3] + '_evalue'] = alignment_values[1]
                        line[sys.argv[3] + '_description'] = Organism_data[0]
                        line[sys.argv[3] + '_position'] = str(alignment_values[3])
                        line[sys.argv[3] + '_function'] = Organism_data[1]
                    else:
                        line[sys.argv[3] + '_id'] = ''
                        line[sys.argv[3] + '_evalue'] = ''
                        line[sys.argv[3] + '_description'] = ''
                        line[sys.argv[3] + '_position'] = ''
                        line[sys.argv[3] + '_function'] = ''
                else:
                    pass



elif sys.argv[2].lower() == 'dtaselect':
    line = __MSParser__.dtaselect(item, organism_list)
else:
    print("What type of file are you searching: cimage or dtaselect?")

# file output
output = open(path + 'Cysteine_annotation.txt', 'w')
output.write(__MSParser__.output(header[0], organism_list, sys.argv[3]))
for peptide in peptides:
    output.write(__MSParser__.output(peptide, organism_list, sys.argv[3]))
    for protein in proteins:
        if protein['index'].strip() == peptide['index'].strip():
            output.write(__MSParser__.output(protein, organism_list, sys.argv[3]))
        else:
            pass
output.close()

