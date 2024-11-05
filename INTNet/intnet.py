"""Main module."""
import argparse
import textwrap
#import intnet_lsaa as lsaa
#import intnet_lsnt as lsnt
#import intnet_ssaa as ssaa
#import intnet_ssnt as ssnt
import sys
parser = argparse.ArgumentParser(
prog='INTNet',
formatter_class=argparse.RawDescriptionHelpFormatter,
description=textwrap.dedent("""\
    IntNet: multi-task multilabel deep neural networks for identification and classification of integrons.
   --------------------------------------------------------------------------------------------------------
    The standlone program is at https:...
    The online service is at https:...

    The input can be long amino acid sequences(full length/contigs), long nucleotide sequences,
    short amino acid reads (30-50aa), short nucleotide reads (100-150nt) in fasta format.
    If your input is short reads you should assign 'intnet-s' model, or if your input is full-length/contigs
    you should assign 'intnet-l' to make the predict.

    USAGE:
        for full-length or contigs
            python intnet.py --input input_path_data --type aa/nt --model intnet-l  --outname output_file_name
        for short reads
            python intnet.py --input input_path_data --type aa/nt --model intnet-s  --outname output_file_name

    general options:
        --input/-i    the test file as input
        --type/-t     molecular type of your test data (aa for amino acid, nt for nucleotide)
        --model/-m    the model you assign to make the prediction (intnet-l for long sequences, intnet-s for short reads)
        --outname/-on  the output file name
    """

),
epilog='Hope you enjoy INTNet journey, any problem please contact scpeiyao@gmail.com')

parser.print_help()

parser.add_argument('-i', '--input', required=True, help='the test data as input')
parser.add_argument('-t', '--type', required=True, choices=['aa', 'nt'], help='molecular type of your input file')
parser.add_argument('-m', '--model', required=True, choices=['intnet-s', 'intnet-l'], help='the model to make the prediction')
parser.add_argument('-on', '--outname', required=True, help='the name of results output')

args = parser.parse_args()

## for AESS_aa -> classifier
if args.type == 'aa' and args.model == 'intnet-s':
    import intnet_ssaa as ssaa
    ssaa.intnet_ssaa(args.input, args.outname)

# for AESS_nt -> classifier
if args.type == 'nt' and args.model == 'intnet-s':
    import  intnet_ssnt as ssnt
    ssnt.intnet_ssnt(args.input, args.outname)

# for AELS_aa -> classifier
if args.type == 'aa' and args.model == 'intnet-l':
    import intnet_lsaa_g3 as lsaa
    lsaa.intnet_lsaa(args.input, args.outname)

# for AELS_nt -> classifier
if args.type == 'nt' and args.model == 'intnet-l':
    import intnet_lsnt_g3 as lsnt
    lsnt.intnet_lsnt(args.input, args.outname)
