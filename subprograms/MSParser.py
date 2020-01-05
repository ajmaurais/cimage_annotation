#! /usr/bin/env/ python

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
    line_dict = ''
    return line_dict

def output(line_dict, organism_list, defined_organism):
    print(len(line_dict))
    line = ''
    line = (line_dict['index'] + '\t')
    line += (line_dict['id'] + '\t')
    line += (line_dict['symbol'] + '\t')
    line += (line_dict['description'] + '\t')
    line += (line_dict['protein_location'] + '\t')
    line += (line_dict['sequence'] + '\t')
    line += (line_dict['mass'] + '\t')
    line += (line_dict['position'] + '\t')
    line += (line_dict['cys_function'] + '\t')
    for organism in organism_list:
        line += (line_dict[organism + '_conserved'] + '\t')
    line += (line_dict[defined_organism + '_id'] + '\t')
    line += (line_dict[defined_organism + '_evalue'] + '\t')
    line += (line_dict[defined_organism+ '_description'] + '\t')
    line += (line_dict[defined_organism+ '_position'] + '\t')
    line += (line_dict[defined_organism+ '_function'] + '\t')


    n = 1
    while n < (len(line_dict) - 19):
        print(line_dict[n])
        line += (line_dict[n] + '\t')
        n += 1
    line = line.rstrip('\t')
    return line
