
import sys
import argparse
import re
from multiprocessing import cpu_count

from .submodules import MSParser, UniProt, Alignments, fasta, parent_parser

PROG_VERSION = 2.1
SEQ_PATH = 'sequences.fasta'
FXN_SEP = '!'
RESIDUE_SEP = '|'

def read_input(args):

    if args.file_type == 'cimage':
        ret = MSParser.Cimage_file()
    elif args.file_type == 'tsv':
        ret = MSParser.Tsv_file(id_col=args.id_col, seq_col=args.seq_col)
    elif args.file_type == 'dtaselect':
        ret = MSParser.Dtaselect()
    else:
        raise RuntimeError('{} is an unknown input file_type'.format(args.file_type))

    ret.read(args.input_file, args.defined_organism)
    return ret


def main():
    parser = argparse.ArgumentParser(prog='cimage_annotation', parents=[parent_parser.PARENT_PARSER],
                                     description='Annotate functional cysteine residues in cimage output.',
                                     epilog='cimage_annotation was written by Dan Bak and Aaron Maurais.\n')

    parser.add_argument('-p', '--parallel', choices=[0, 1], type=int, default=1,
                        help='Choose whether internet queries and protein alignments should be performed in parallel.'
                        ' Parallel processing is performed on up to the number of logical cores on your system. '
                        '1 is the default.')

    parser.add_argument('-t', '--nThread', type=int, default=None,
                        help='Chose how many threads to use for parllel processing.'
                        'This option overrides the --parallel option.')

    parser.add_argument('--debug', choices=['none', 'pdb', 'pudb'], default='none',
                        help='Start the main method in the selected debugger.')

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
        _nThread=1


    # Open input file
    input_file = read_input(args)

    if len(input_file) == 0:
        sys.stderr.write('ERROR: No peptides found in {}!\n\tExiting...\n'.format(args.input_file))
        return -1

    sys.stdout.write('\nRetreiving protein Uniprot records...\n')
    record_dict = UniProt.get_uniprot_records(input_file.unique_ids, _nThread, verbose=args.verbose,
            show_bar = not(args.verbose and args.parallel == 0))

    sequences = dict()
    seq_written = False
    for i, p in input_file.iterpeptides():
        # Get and Parse Uniprot entry for protein
        try:
            seq_temp = re.match(r'.?\.?([A-z\*]+)\.?/?', p[args.seq_col]).group(1)
        except AttributeError as e:
            sys.stdout.write('Error parsing sequence: {}'.format(p[args.seq_col]))
            continue

        UniProt_data = UniProt.ExPasy(seq_temp, record_dict[p[args.id_col]],
                                      all_features=args.all_features, res_sep=RESIDUE_SEP, fxn_sep=FXN_SEP)

        input_file.set_peptide_value(i, 'position', UniProt_data[1])   # cysteine position
        input_file.set_peptide_value(i, 'res_function', UniProt_data[2]) # cysteine function (if known)
        input_file.set_peptide_value(i, 'domains', UniProt_data[3]) # Domain at position (if known)
        input_file.set_peptide_value(i, 'protein_location', UniProt_data[5]) # protein subcellular localization (if known)

        if args.write_seq:
            if p[args.id_col] not in sequences: #only write sequence if it is not currently in file
                fasta.write_fasta_entry(SEQ_PATH,
                                        p[args.id_col],
                                        UniProt_data[3],
                                        description=p['description'],
                                        append=seq_written)
                seq_written = True
        sequences[p[args.id_col]] = ('' if args.description_col not in p else p[args.description_col], UniProt_data[3])

    # Replace alignment files with empty string so they won't be continuously appended to.
    if args.write_alignment_data and args.align:
        # blast protein sequence against each fasta database and parse those alignments to determine cysteine conservation
        sys.stdout.write('\nAlligning protein sequences to determine cysteine conservation...\n')
        for organism in Alignments.organism_list:
            align_fname = '{}_alignments.{}'.format(organism, args.align_format)
            sys.stdout.write('\tCreating {}...'.format(align_fname))
            with open(align_fname, 'w') as outF:
                outF.write('')
            sys.stdout.write('Done!\n')

    if args.align:
        alignment_data = Alignments.align_all(input_file.unique_ids, sequences, args.database_dir, Alignments.organism_list,
                                              nThread=_nThread, verbose=args.verbose,
                                              show_bar=not(args.verbose and args.parallel == 0))
        
        for o in Alignments.organism_list:
            input_file.add_column('{}_conserved'.format(o))

        if args.defined_organism != 'none':
            
            for o in MSParser.ALLIGNMENT_COLUMNS:
                input_file.add_column('{}_{}'.format(args.defined_organism, o))

            org_ids = set()
            for a in alignment_data.values():
                org_ids.add(a[args.defined_organism].get_best_id())

            sys.stdout.write('\nRetreiving Uniprot records for {} alignments...\n'.format(args.defined_organism))
            org_record_dict = UniProt.get_uniprot_records(org_ids, _nThread, verbose=args.verbose,
                                                          show_bar=not(args.verbose and args.parallel==0))

        seen=set() # Keep track of seen protein IDs
        for i, p in input_file.iterpeptides():
            for organism in Alignments.organism_list:
                if p[args.id_col] in seen and args.write_alignment_data:
                    alignment_data[p[args.id_col]][organism].write('{}_alignments.{}'.format(organism, args.align_format),
                                                            file_format=args.align_format, mode='a')
                seen.add(p[args.id_col])

                evalue = alignment_data[p[args.id_col]][organism].get_best_evalue()
                conserved_temp = list()
                for pos in p['position'].split(RESIDUE_SEP):
                    cp_temp = '--'
                    if pos in ('BAD_ID', 'RESIDUE_NOT_FOUND'):
                        cp_temp = 'Error'
                    else:
                        assert(pos.isdigit())
                        if evalue is None:
                            cp_temp == '--'
                        elif evalue <= args.evalue_co:
                            cp_temp = 'Yes' if alignment_data[p[args.id_col]][organism].conserved_at_position(int(pos)) else 'No'
                    conserved_temp.append(cp_temp)

                input_file.set_peptide_value(i, '{}_conserved'.format(str(organism)), RESIDUE_SEP.join(conserved_temp))

                # for comparative organism analyze Uniprot entry of best blast hit
                if organism == args.defined_organism.lower():
                    org_dict_temp = {x: '' for x in MSParser.ALLIGNMENT_COLUMNS}

                    id_temp = alignment_data[p[args.id_col]][organism].get_best_id()
                    org_dict_temp['id'] = id_temp
                    org_dict_temp['description'] = alignment_data[p[args.id_col]][organism].get_best_description()
                    if id_temp != '':
                        org_dict_temp['evalue'] = alignment_data[p[args.id_col]][organism].get_best_evalue()
                        if org_dict_temp['evalue'] <= args.evalue_co:

                            positions_temp = list()
                            functions_temp = list()
                            for pos in p['position'].split(RESIDUE_SEP):
                                homolog_position = None if not pos.isdigit() else alignment_data[p[args.id_col]][organism].alignment_at_position(int(pos))[1]
                                if homolog_position is None:
                                    homolog_position = 0
                                positions_temp.append(str(homolog_position))
                                if org_record_dict[id_temp] is None:
                                    functions_temp.append('')
                                else:
                                    functions_temp.append(UniProt.res_features(org_record_dict[id_temp],
                                                                               homolog_position - 1,
                                                                               all_features=args.all_features))

                        org_dict_temp['position'] = RESIDUE_SEP.join(positions_temp)
                        if ''.join(functions_temp):
                            if len(functions_temp) == 1:
                                org_dict_temp['function'] = functions_temp[0]
                            else:
                                org_dict_temp['function'] = FXN_SEP.join(['{}:{}'.format(p, s) for p, s in zip(positions_temp,
                                                                                                              functions_temp)])

                    # add alignment data to peptides
                    for k, v in org_dict_temp.items():
                        key_temp = '{}_{}'.format(args.defined_organism, k)
                        input_file.set_peptide_value(i, key_temp, v)

    # file output
    input_file.write(args.ofname)
    sys.stdout.write('\nResults written to {}\n\n'.format(args.ofname))


if __name__ == '__main__':
    main()

