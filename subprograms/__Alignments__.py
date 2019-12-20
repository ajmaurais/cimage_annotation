#! usr/bin/env/ python

import sys

def align_write(path, organism):
    f = open(path + organism + '_alignments.txt', 'a')
    g = open(path + 'alignment.txt', 'r')
    align = g.read()
    f.write(align + '\n\n\n\n\n\n')
    f.close()
    g.close()

def Align_file_parse(position, path):
    f = open(path + 'alignment.txt', 'r')
    align = f.read()
    f.close()

    conserved = ''
    id = ''
    evalue = ''
    homolog_position = ''
    seq_list = []
    if len(align.split('\n\n')) < 9:
        conserved = 'Error'
    else:
        align_header = align.split('\n\n')[8]

        if align_header.strip() == '***** No hits found *****':
            conserved = '--'
        else:
            aligned_proteins = align_header.split('\n')

            id = aligned_proteins[0].split('|')[1]
            evalue = aligned_proteins[0].strip().split(' ')[-1]
        
            if float(evalue) > float('1e-5'):
                conserved = '--'
            else:
            
                align_data = align.split('>')[1].split('\n\n\n')[0].split('\n\n')
            
                for item in align_data:
                    if item[:5] == 'Query':
                        query = item.split('\n')[0].split()
                        subject = item.split('\n')[2].split()

                        if int(position) >= int(query[1]) and int(position) <= int(query[3]):
                            m = int(query[1])
                            n = int(subject[1])
                            i = 0
                            while i < len(query[2]):
                                if query[2][i] == '-':
                                    a = m
                                else:
                                    a = m
                                    m += 1
                                if subject[2][i] == '-':
                                    b = n
                                else:
                                    b = n
                                    n += 1
                                seq_list.append(str(a) + '\t' + query[2][i] + '\t' + str(b) + '\t' + subject[2][i])

                                i += 1
                        else:
                            pass
                    else:
                        pass

            if seq_list == []:
                conserved = '--'
            else:
                for item in seq_list:
                    if item.split('\t')[0] == str(position) and item.split('\t')[1] == item.split('\t')[3]:
                        conserved = 'Yes'
                        homolog_position = item.split('\t')[2]
                    else:
                        pass
            if homolog_position == '' and conserved == '':
                conserved = 'No'

    return id, evalue, conserved, homolog_position
