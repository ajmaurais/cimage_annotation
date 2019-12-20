#! /usr/bin/env/ python

import sys
import os
import os.path
import urllib,urllib2
import httplib
sys.path.insert(0, 'biopython-1.65')

from Bio import ExPASy, SwissProt, SeqIO

features_list =['CA_BIND', 'ZN_FING', 'DNA_BIND', 'NP_BIND',
                 'ACT_SITE', 'METAL', 'BINDING', 'SITE',
                 'NON_STD', 'MOD_RES', 'LIPID', 'CARBOHYD',
                 'DISULFID', 'CROSSLINK', 'VARIANT', 'MUTAGEN',
                 'UNSURE', 'CONFLICT', 'REGION']

def protein_location(record):
    pro_location = ''
    for item in record.comments:

        if item.split(':')[0] == 'SUBCELLULAR LOCATION':
            if item.split(':')[1].strip().split(' ')[0] == 'Isoform':
                pro_location += item.split(':')[2]
            else:
                pro_location += item.split(':')[1]
        else:
            pass
    return str(pro_location)

    
def cys_position(record, sequence, location):
    peptide_location = str(record.sequence).find(sequence)
    if str(peptide_location) == str(-1):
        cys_location = 0
    else:
        cys_location = peptide_location + location

    return str(cys_location)


def cys_function(record, position):
    cys_function = ''

    for feature in record.features:

        if feature[0].upper() == 'DISULFID' and str(feature[1])[0] != '?' and str(feature[2])[0] != '?' and int(position) >= int(feature[1]) and int(position) <= int(feature[2]):
            cys_function += (feature[0] + '--' + feature[3] + ' || ')
        elif feature[0].upper() in features_list and str(feature[1])[0] != '?' and str(feature[2])[0] != '?' and (int(feature[2]) - int(feature[1])) <= 10 and int(position) >= int(feature[1]) and int(position) <= int(feature[2]):
            cys_function += (feature[0] + '--' + feature[3] + ' || ')
    return cys_function

def ExPasy(id, sequence, location):
    position = ''
    function = ''
    organism = ''
    full_sequence = ''
    pro_location = ''
    try:
        handle = ExPASy.get_sprot_raw(id)
        try:
            record = SwissProt.read(handle)
        except ValueError:
            record = None
    except (urllib2.HTTPError, httplib.BadStatusLine, urllib2.URLError), err:
        print err
	record = None
        
    if record is not None:
        organism = record.organism
        pro_location = protein_location(record)
        position = cys_position(record, sequence, location)
        function = cys_function(record, position)
        full_sequence = record.sequence
    else:
        position = 'Bad ID'
    return organism, position, function, full_sequence, pro_location

def ExPasy_alt(id, position):
    protein = ''
    function = ''

    try:
        handle = ExPASy.get_sprot_raw(id)
        try:
            record = SwissProt.read(handle)
        except ValueError:
            record = None
    except (urllib2.HTTPError, httplib.BadStatusLine, urllib2.URLError), err:
        print err
	record = None

    if record is not None:
        protein = record.description
        function = cys_function(record, position)

    else:
        protein = 'Bad ID'

    return protein, function
