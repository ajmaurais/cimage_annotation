
import os
import os.path
import sys
import argparse
import re
from multiprocessing import cpu_count

from .submodules import MSParser, UniProt, Blast, Alignments, Fasta, parent_parser

PROG_VERSION=2.0
SEQ_PATH='sequences.fasta'

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
    parser = argparse.ArgumentParser(prog = 'cimage_annotation', parents=[parent_parser.PARENT_PARSER],
                                     description = 'Annotate functional cysteine residues in cimage output.')

    parser.add_argument('-p', '--parallel', choices=[0,1], type=int, default=1,
                        help='Choose whether internet queries and protein alignments should be performed in parallel.'
                        ' Parallel processing is performed on up to the number of logical cores on your system. '
                        '1 is the default.')

    parser.add_argument('-t', '--nThread', type=int, default=None,
                        help='Chose how many threads to use for parllel processing.'
                        'This option overrides the --parallel option.')

    parser.add_argument('--debug', choices=['none', 'pdb', 'pudb'], default='none',
            help='Start the main method in the selected debuger.')

    args = parser.parse_args()

    if args.debug != 'none':
        assert(args.debug in ['pdb', 'pudb'])
        if args.debug == 'pdb':
            import pdb as db
        elif args.debug == 'pudb':
            try:
                import pudb as db
            except ModuleNotFoundError as e:
                sys.stderr.write('pudb is not installed.')
                return -1
        db.set_trace()

    sys.stdout.write('\ncimage_annotation v{}\n'.format(PROG_VERSION))

    # Manually check args.
    if args.align and args.database_dir is None:
        sys.stderr.write('--database_dir must be specified when --align 1 is set\n')
        return -1
    _nThread = args.nThread
    if args.parallel and args.nThread is None:
        _nThread = cpu_count()
    elif not args.parallel and args.nThread is None:
        _nThread = 1


    # Open cimage or dtaselect file
    file_input = parse_input(args.input_file, args.file_type, args.defined_organism)

    header = []
    cysteines = []
    peptides = []
    index = ''
    uniuque_ids = set()
    sequences = dict()

    # Analysis for a cimage file creates lists (above) of dictionaries
    for i, line in enumerate(file_input):

        if line['index'].strip() == 'index':  # create header line(s)
            header.append(line)

        elif line['index'].strip() != '':  # create peptide lines (no protein information - includes overall ratio)
            index = line['index'].strip()
            cysteines.append(line)

        else:
            file_input[i]['index'] = index # create protein lines
            peptides.append(line)
            uniuque_ids.add(line['id'])

    if len(peptides) == 0:
        sys.stderr.write('ERROR: No peptides found in {}!\n\tExiting...\n'.format(args.input_file))
        return -1

    sys.stdout.write('\nRetreiving protein Uniprot records...\n')
    record_dict = UniProt.get_uniprot_records(uniuque_ids, _nThread, verbose=args.verbose,
            show_bar=not(args.verbose and args.parallel==0))

    seq_written=False
    for i, p in enumerate(peptides):
        # Get and Parse Uniprot entry for protein
        try:
            seq_temp = re.match('.?\.?([A-z\*]+)\.?/?', p['sequence']).group(1)
        except AttributeError as e:
            sys.stdout.write('Error parsing sequence: {}'.format(p['sequence']))
            continue

        UniProt_data = UniProt.ExPasy(seq_temp, record_dict[p['id']],
                                      res_sep=args.residue_sep, fxn_sep=args.fxn_sep)
        peptides[i]['position'] = UniProt_data[1]   # cysteine position
        peptides[i]['cys_function'] = UniProt_data[2] # cysteine function (if known)
        peptides[i]['protein_location'] = UniProt_data[4] # protein subcellular localization (if known)

        if args.write_seq:
            if p['id'] not in sequences: #only write sequence if it is not currently in file.
                Fasta.write_fasta_entry(SEQ_PATH,
                                        p['id'],
                                        UniProt_data[3],
                                        description=p['description'],
                                        append=seq_written)
                seq_written = True
        sequences[peptides[i]['id']] = (p['description'], UniProt_data[3])

    # Replace alignment files with empty string so they won't be continuously appended to.
    if args.write_alignment_data and args.align:
        # blast protein sequence against each fasta database and parse those alignments to determine cysteine conservation
        sys.stdout.write('\nAlligning protein sequences to determine cysteine conservation...\n')
        for organism in organism_list:
            align_fname = '{}_alignments.{}'.format(organism, args.align_format)
            sys.stdout.write('\tCreating {}...'.format(align_fname))
            with open(align_fname, 'w') as outF:
                outF.write('')
            sys.stdout.write('Done!\n')

    if args.align:
        alignment_data = Alignments.align_all(peptides, sequences, args.database_dir, organism_list,
                                              nThread=_nThread, verbose=args.verbose,
                                              show_bar=not(args.verbose and args.parallel==0))

        if args.defined_organism != 'none':
            org_ids = set()
            for a in alignment_data.values():
                org_ids.add(a[args.defined_organism].get_best_id())

            sys.stdout.write('\nRetreiving Uniprot records for {} alignments...\n'.format(args.defined_organism))
            org_record_dict = UniProt.get_uniprot_records(org_ids, _nThread, verbose=args.verbose,
                    show_bar=not(args.verbose and args.parallel==0))

        seen=set() # Keep track of seen protein IDs
        for i, p in enumerate(peptides):
            for organism in organism_list:
                if p['id'] in seen and args.write_alignment_data:
                    alignment_data[p['id']][organism].write('{}_alignments.{}'.format(organism, args.align_format),
                                                            file_format = args.align_format, mode = 'a')
                seen.add(p['id'])

                evalue = alignment_data[p['id']][organism].get_best_evalue()
                conserved_temp = list()
                for pos in p['position'].split(args.residue_sep):
                    cp_temp = '--'
                    if pos == 'BAD_ID':
                        cp_temp = 'Error'
                    else:
                        assert(pos.isdigit())
                        assert(evalue is not None)
                        if evalue <= args.evalue_co:
                            cp_temp = 'Yes' if alignment_data[p['id']][organism].conserved_at_position(int(pos)) else 'No'
                    conserved_temp.append(cp_temp)

                peptides[i][str(organism) + '_conserved'] = args.residue_sep.join(conserved_temp)

                # for comparative organism analyze Uniprot entry of best blast hit
                if organism == args.defined_organism.lower():
                    org_dict_temp = {x: list() for x in ['id', 'evalue', 'description', 'position', 'function']}

                    id_temp = alignment_data[p['id']][organism].get_best_id()
                    org_dict_temp['id'].append(id_temp)
                    org_dict_temp['description'].append(alignment_data[p['id']][organism].get_best_description())
                    if id_temp != '':
                        org_dict_temp['evalue'].append(alignment_data[p['id']][organism].get_best_evalue())
                        if org_dict_temp['evalue'][0] <= args.evalue_co:
                            for pos in p['position'].split(args.residue_sep):
                                homolog_position = alignment_data[p['id']][organism].alignment_at_position(int(pos))[1]
                                org_dict_temp['position'].append(homolog_position)

                                org_dict_temp['function'].append(UniProt.cys_function(org_record_dict[id_temp], homolog_position - 1))

                    # concatenate alignment data
                    org_dict_temp = {k: args.residue_sep.join([str(x) for x in v]) for k, v in org_dict_temp.items()}

                    # add alignment data to peptides
                    for k, v in org_dict_temp.items():
                        key_temp = '{}_{}'.format(args.defined_organism, k)
                        peptides[i][key_temp] = v

    # file output
    with open(args.ofname, 'w') as outF:
        outF.write(MSParser.output(header[0], organism_list, args.defined_organism))
        for cysteine in cysteines:
            outF.write(MSParser.output(cysteine, organism_list, args.defined_organism))
            for peptide in peptides:
                if peptide['index'].strip() == cysteine['index'].strip():
                    outF.write(MSParser.output(peptide, organism_list, args.defined_organism))
    sys.stdout.write('\nResults written to {}\n\n'.format(args.ofname))


if __name__ == '__main__':
    main()

