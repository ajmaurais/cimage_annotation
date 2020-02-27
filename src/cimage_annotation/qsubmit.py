
import sys
import os
import argparse
import subprocess
import os.path
from math import ceil

from .submodules import parent_parser

BLAST_PBS_VERSION = 'blast'
PBS_MODULE_LOAD_COMMAND = 'module load'
CIMAGE_ANNOTATION_EXE = 'cimage_annotation'

def makePBS(mem, ppn, walltime, wd, cimage_annotation_args):
    pbsName = '{}/cimage_annotation.pbs'.format(wd)
    _flags = ' '.join(['--{} {}'.format(k,v) for k, v in cimage_annotation_args.items() if v is not None and k != 'input_file'])

    sys.stdout.write('Writing {}...'.format(pbsName))
    with open(pbsName, 'w') as outF:
        outF.write("#!/bin/tcsh\n")
        outF.write('#PBS -l mem={}gb,nodes=1:ppn={},walltime={}\n\n'.format(mem, ppn, walltime))
        outF.write('{} {}\n\n'.format(PBS_MODULE_LOAD_COMMAND, BLAST_PBS_VERSION))
        outF.write('cd {}\n'.format(wd))
        outF.write('{} {} {} > stdout.txt\n'.format(CIMAGE_ANNOTATION_EXE, _flags, cimage_annotation_args['input_file']))

    sys.stdout.write('Done!\n')
    return pbsName


def main():

    parser = argparse.ArgumentParser(prog='qsub_cimage_annotation', parents=[parent_parser.PARENT_PARSER],
                                     description='Submit cimage_annotation job to the queue.')

    parser.add_argument('-m', '--mem', default=4, type=int,
                        help='Amount of memory to allocate per PBS job in gb. Default is 4.')

    parser.add_argument('-p', '--ppn', default=8, type=int,
                        help='Number of processors to allocate per PBS job. Default is 8.')

    parser.add_argument('-t', '--walltime', default='12:00:00',
                        help='Walltime per job in the format hh:mm:ss. Default is 12:00:00.')

    parser.add_argument('-g', '--go', action='store_true', default=False,
                        help='Should job be submitted? If this flag is not supplied, program will be a dry run. '
                             '.pbs file will written but job will not be submitted.')

    args = parser.parse_args()
    parent_args = parent_parser.PARENT_PARSER.parse_known_args()[0]

    # calc nThread
    if args.parallel and args.nThread is None:
        _nThread = args.ppn * 2
    elif not args.parallel and args.nThread is None:
        _nThread=1
    else:
        _nThread = args.nThread
        args.ppn = ceil(args.nThread / 2) * 2

    # Manually check args
    n_arg_errors = 0
    error_message = '\nThere were errors in the options you specified...'
    if not os.path.isfile(args.input_file):
        error_message += '\n\tSpecified input file: {} does not exist on path:\n\t\t{}\n'.format(args.input_file,
                                                                                                 os.path.abspath(args.input_file))
        n_arg_errors += 1
    if args.align and args.database_dir is None:
        error_message += '\n\t--database_dir must be specified when --align 1 is set\n'
        n_arg_errors += 1
    if n_arg_errors > 0:
        sys.stderr.write(error_message)
        return -1

    sys.stdout.write('\nRequested job with {} processor and {}gb of memory...\n'.format(args.ppn, args.mem))
    #get wd
    wd = os.path.dirname(os.path.abspath(args.input_file))

    cimage_annotation_args = {arg: getattr(args, arg) for arg in vars(parent_args)}
    cimage_annotation_args['nThread'] = _nThread
    cimage_annotation_args['align'] = '' if args.align else None
    cimage_annotation_args['verbose'] = '' if args.verbose else None
    cimage_annotation_args['write_seq'] = '' if args.write_seq else None
    cimage_annotation_args['all_features'] = '' if args.all_features else None

    pbsName = makePBS(args.mem, args.ppn, args.walltime, wd, cimage_annotation_args)
    command = 'qsub {}'.format(pbsName)
    if args.verbose:
        sys.stdout.write('{}\n'.format(command))
    if args.go:
        proc = subprocess.Popen([command], cwd=wd, shell=True)
        proc.wait()

if __name__ == '__main__':
    main()

