#!/usr/bin/env python3

"""
This script update vcf file to add both HP haplotag and PS phasing block info fields, It takes as input vcf file, hp, ps.
"""
import argparse
import sys, os
from operator import itemgetter
from collections import Counter

# Python program to print
# green text with red background

from colorama import init
from termcolor import colored

init()




def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Produce phasing report',
                                     add_help=True, )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
    # parser.add_argument('input', help='Input file', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    # parser.add_argument('output', help='Output file', nargs="?", type=argparse.FileType('w'), default=sys.stdout)

    parser.add_argument('input', nargs='?', help="Structural variant vcf file",
                             type=argparse.FileType('r'),
                             default=sys.stdin)
    parser.add_argument('hp', nargs='?', help="tab delimeted read\thp\tps file",
                             type=argparse.FileType('r'))
    parser.add_argument('output', nargs='?', help="Output file, PS and HP will be added.",
                                 type=argparse.FileType('w+'),
                                 default=sys.stdout)

    parser.set_defaults(func=update_vcf)

    # if not argument print help.
    if len(sys.argv) == 1 and  sys.stdin.isatty():  # sys.stdin.isatty() returns false if there's something in stdin
         parser.print_help(sys.stderr)
         sys.exit(1)

    args = parser.parse_args()


    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

def update_vcf(args):
    # check if the input from stdin
    if not sys.stdin.isatty(): # there is nothing in the stdin
        if args.input.name.endswith("gz"):
            import gzip
            myfile = gzip.open(args.input.name, 'rt') # t is not a must normally it is default.
        else:
            myfile = args.input
    else:
        myfile = args.input

    # read the SV file as dictionary
    hp_dic = {}
    with args.hp as hp_in:
        for line in hp_in:
            id, hp, ps = line.split()
            hp_dic[id] = [hp.rsplit(":", 1)[-1], ps.rsplit(":", 1)[-1]]

    with myfile as data_in, args.output as data_out:
        for line in data_in:
            reads = []
            if line.startswith('##'):
                data_out.write(line)
            elif line.startswith("#"):
                # data_out.write("##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Haplotype identifier\">\n")
                data_out.write("##INFO=<ID=CONFLICT,Number=1,Type=Integer,Description=\"The Phase is conflict or not\">\n")
                data_out.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n")
                data_out.write(line)
            else:
                line_split = line.split()
                if line_split[-1].split(":", 1)[0] == "1/1" or line_split[-1].split(":", 1)[0] == "0/0":  # no gt to phase
                    data_out.write("{}\n".format("\t".join(line_split)))
                else:
                    reads = line_split[7].split(";")[10].split(",")
                    reads[0] = reads[0].split("=")[-1]
                    myvalues = list(map(hp_dic.get, reads))  # list of lists first element id hp second is ps or Nonn on case there are no reads with hp and ps to support this sv

                    # If not all None
                    if any(myvalues): # any value is not none
                        ps_dict = categorize_ps(myvalues)
                        if 0 in list(ps_dict.values()): # means that the hp is conflicting do not update anything and add flag that is is conflicting.
                            line_split[7] = "{info};CONFLICT={conflict}".format(info=line_split[7], conflict=1)
                            line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                            line_split[-1] = "{}:{}".format(line_split[-1], ",".join(ps_dict.keys()))
                            data_out.write("{}\n".format("\t".join(line_split)))
                        else: # update the gt field and ps to sv
                            line_split[7] = "{info};CONFLICT={conflict}".format(info=line_split[7], conflict=0)
                            line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                            line_split[-1] = line_split[-1].replace("/", "|")
                            line_split[-1] = "{}:{}".format(line_split[-1], ",".join(ps_dict.keys()))
                            data_out.write("{}\n".format("\t".join(line_split)))
                    else: # all are none
                        line_split[7] = "{info};CONFLICT=2".format(info=line_split[7])


def categorize_ps(myvalues):
    myvalues = [i for i in myvalues if i is not None] # remove None
    ps_dict = {}
    for i in myvalues:
        ps = (i[1])
        hp = int(i[0])
        if ps in ps_dict:
            if int(hp) == 1:
                if ps_dict[ps] < 0 :
                    ps_dict[ps] = ps_dict[ps] - 1
                else: #conflict
                    ps_dict[ps] = 0

            else: # means that it is haplotype 2 hp=2
                if ps_dict[ps] > 0:
                    ps_dict[ps] = ps_dict[ps] + 1
                else: #conflict
                    ps_dict[ps] = 0
        else:
            if int(hp) == 1:
                ps_dict[ps] = -1
            else:
                ps_dict[ps] = 1
    return ps_dict


def main():
    args = get_args()



if __name__ == "__main__":
    main()
