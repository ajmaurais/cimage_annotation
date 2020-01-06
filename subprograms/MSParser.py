#! /usr/bin/env/ python

PRINT_COLS=['index', 'id', 'symbol', 'description', 'protein_location', 'sequence', 'mass', 'position', 'cys_function']
ALLIGNMENT_COLUMNS = ['id', 'evalue', 'description', 'position', 'function']

def cimage(line, organism_list, defined_organism):
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
        line_dict['position'] = 'cysteine position'
        line_dict['cys_function'] = 'cysteine function'
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
        line_dict['cys_function'] = ''
        line_dict[defined_organism + '_id'] = ''
        line_dict[defined_organism + '_evalue'] = ''
        line_dict[defined_organism + '_description'] = ''
        line_dict[defined_organism + '_position'] = ''
        line_dict[defined_organism + '_function'] = ''
        for organism in organism_list:
            line_dict[organism + '_conserved'] = ''

    return line_dict

def dtaselect(line, organism_list):
    raise NotImplementedError()

    line_dict = ''
    return line_dict

def output(line_dict, organism_list, defined_organism):

    line = '\t'.join([line_dict[x] for x in PRINT_COLS])

    for organism in organism_list:
        line += (line_dict[organism + '_conserved'] + '\t')
    line += '\t'.join([line_dict['{}_{}'.format(defined_organism, x)] for x in ALLIGNMENT_COLUMNS])

    n = 1
    while n < (len(line_dict) - 19):
        line += (line_dict[n] + '\t')
        n += 1
    line = line.rstrip('\t')
    return line
