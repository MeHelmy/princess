#!/usr/bin/env python3

"""
This script update Methylation file to add both HP haplotag and PS phasing block, It takes as input meth file, hp, ps.
"""
import argparse
import sys, os
from operator import itemgetter
from collections import Counter

# Python program to print
# green text with red background
#
# from colorama import init
# from termcolor import colored
#
# init()




def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Produce phasing report for Methylation',
                                     add_help=True, )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
    # parser.add_argument('input', help='Input file', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    # parser.add_argument('output', help='Output file', nargs="?", type=argparse.FileType('w'), default=sys.stdout)

    parser.add_argument('input', nargs='?', help="Methylation file",
                             type=argparse.FileType('r'),
                             default=sys.stdin)
    parser.add_argument('hp', nargs='?', help="tab delimeted read\thp\tps file",
                             type=argparse.FileType('r'))
    parser.add_argument('output', nargs='?', help="Output file, PS and HP will be added.",
                                 type=argparse.FileType('w+'),
                                 default=sys.stdout)

    parser.set_defaults(func=update_meth)

    # if not argument print help.
    if len(sys.argv) == 1 and  sys.stdin.isatty():  # sys.stdin.isatty() returns false if there's something in stdin
         parser.print_help(sys.stderr)
         sys.exit(1)

    args = parser.parse_args()


    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

def update_meth(args):
    # check if the input from stdin
    if not sys.stdin.isatty(): # there is nothing in the stdin
        if args.input.name.endswith("gz"):
            import gzip
            myfile = gzip.open(args.input.name, 'rt') # t is not a must normally it is default.
        else:
            myfile = args.input
    else:
        myfile = args.input

    # read the Haplotype file as dictionary
    hp_dic = {}
    with args.hp as hp_in:
        for line in hp_in:
            id, hp, ps = line.split()
            hp_dic[id] = [hp.rsplit(":", 1)[-1], ps.rsplit(":", 1)[-1]] # read hp, ps


    with myfile as data_in, args.output as data_out:
        first = True
        n = 0
        for line in data_in:
            n+=1
            if first:
                first = False
                data_out.write(line.strip()+"\tHP\tPS\n")
                continue
            line_split = line.split()
            read = line_split[4]
            hp, ps = hp_dic.get(read, ['.', '.']) # In case f the read have not been haplotyped.
            data_out.write("{}\t{}\t{}\n".format(line.strip(), hp, ps))

def main():
    args = get_args()



if __name__ == "__main__":
    main()
