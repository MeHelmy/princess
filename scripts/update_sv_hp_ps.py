
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
#
# from colorama import init
# from termcolor import colored
#
# init()




def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Phase SVs Using haplotyped reads in tab format',
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
    parser.add_argument('-c', '--conflict', dest='ignore_conflict', metavar='Max Conflict Reads', type=int,  help='Minumum number of conflict reads to ignore', default=0)

    parser.set_defaults(func=update_vcf)

    # if no argument print help.
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

    # read the Haplotyped reads file as dictionary
    hp_dic = {}
    with args.hp as hp_in:
        for line in hp_in:
            id, hp, ps = line.split()
            hp_dic[id] = [hp.rsplit(":", 1)[-1], ps.rsplit(":", 1)[-1]] # read -> [hp, ps]


    with myfile as data_in, args.output as data_out:
        for line in data_in:
            reads = []
            if line.startswith('##'):
                data_out.write(line)
            elif line.startswith("#"):
                # data_out.write("##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Haplotype identifier\">\n")
                data_out.write("##INFO=<ID=CONFLICT,Number=.,Type=Integer,Description=\"The Phase is conflict [1], not [0] or no data [2]\">\n")
                data_out.write("##FORMAT=<ID=PS,Number=.,Type=Integer,Description=\"Phase set identifier\">\n")
                data_out.write(line)
            else:
                line_split = line.split()
                if line_split[-1].split(":", 1)[0] == "1/1" or line_split[-1].split(":", 1)[0] == "0/0" or line_split[-1].split(":", 1)[0] == "./.":  # no gt to phase
                    data_out.write("{}\n".format("\t".join(line_split)))
                elif line_split[-1].split(":", 1)[0] == "0/1" or line_split[-1].split(":", 1)[0] == "1/0":
                    reads = [i for i in line_split[7].split(";")  if i.startswith("RNAMES")][0].split("=",1)[-1].split(",")
                    svtype = [i for i in line_split[7].split(";")  if i.startswith("SVTYPE")][0].split("=",1)[-1].split(",")
                    svlen = [i for i in line_split[7].split(";")  if i.startswith("SVLEN")][0].split("=",1)[-1].split(",")
                    #reads = line_split[7].split(";")[10].split(",") #info field -> reads
                    #reads[0] = reads[0].split("=")[-1]
                    myvalues = list(map(hp_dic.get, reads))  # list of lists first element id hp second is ps or None on case there are no reads with hp and ps to support this sv
                    id = line_split[2]
                    # If any value not None
                    # print(f'{line_split[1]}\t{id}\t{svtype[0]}\t{svlen[0]}\t{myvalues}')
                    if any(myvalues): # any value is not none
                        # print(f'{id}\t{svtype[0]}\t{svlen[0]}\t{myvalues}')
                        ps_dict = categorize_ps_up(myvalues, args.ignore_conflict)
                        if 0 in list(ps_dict.values()): # means that the hp is conflicting do not update anything and add flag that is is conflicting.
                            line_split[7] = "{info};CONFLICT={conflict}".format(info=line_split[7], conflict=1)
                            line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                            line_split[-1] = "{}:{}".format(line_split[-1], ",".join(ps_dict.keys()))
                            data_out.write("{}\n".format("\t".join(line_split)))
                        else: # update the gt field and ps to sv
                            line_split[7] = "{info};CONFLICT={conflict}".format(info=line_split[7], conflict=0)
                            line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                            # if values are negative then it is hp=1 1|0 else it is hp2 0|1
                            # line_split[-1] = line_split[-1].replace("/", "|")
                            hp_new_value = line_split[-1].split(':')
                            try:
                                if list(ps_dict.values())[0] < 1: # haplotype 1
                                    hp_new_value[0] = "1|0"
                                else:
                                    hp_new_value[0] = "0|1"
                            except Exception as e:
                                print(e)


                            hp_new_value = ":".join(hp_new_value)
                            line_split[-1] = "{}:{}".format(hp_new_value, ",".join(ps_dict.keys()))
                            data_out.write("{}\n".format("\t".join(line_split)))
                    else: # all are none
                        line_split[7] = "{info};CONFLICT=2".format(info=line_split[7])
                        line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                        line_split[-1] = "{}:{}".format(line_split[-1], ".")
                        data_out.write("{}\n".format("\t".join(line_split)))


# Test case [['1', '23200'], ['2', '23200'], ['2', '23200'], ['1', '23200'], ['2', '23200'], ['2', '23200'], ['1', '23200'], ['2', '23200'], ['1', '23200'], ['2', '23200'], ['1', '23200'], ['1', '23200'], ['2', '23200'], ['2', '23200'], ['2', '23200']]
#
# [['1', '13164067'], ['1', '13164067'], ['1', '13164067'], ['1', '13164067'], ['1', '13164067'], ['2', '12948612'], ['1', '13164067'], ['1', '13164067'], ['2', '12948612'], ['1', '13164067'], ['2', '12948612']]
def categorize_ps(myvalues):
    myvalues = [i for i in myvalues if i is not None] # remove None
    ps_dict = {}
    for i in myvalues:
        hp = int(i[0])
        ps = i[1]
        if ps in ps_dict:
            if hp == 1:
                if ps_dict[ps] < 0:
                    ps_dict[ps] = ps_dict[ps] - 1
                else: #conflict
                    ps_dict[ps] = 0

            else: # means that it is haplotype 2 hp=2
                if ps_dict[ps] > 0:
                    ps_dict[ps] = ps_dict[ps] + 1
                else: #conflict
                    ps_dict[ps] = 0
        else:
            if hp == 1:
                ps_dict[ps] = -1
            else:
                ps_dict[ps] = 1
    return ps_dict


def most_frequent(List):
	return max(set(List), key = List.count)


def categorize_ps_conflict(myvalues, max_conflict):
    ps_dict = {}
    myvalues = [i for i in myvalues if i is not None] # remove None
    hp = [i[0] for i in myvalues]
    hp_count = {i: hp.count(i) for i in hp} # i.e {'1': 3, '2': 2} or {'1': 3}
    # TODO: chek if we have two hap with differnt phase block they should not be counted as conflict
    if len(hp_count) > 1: # they are conflicting
        if min(hp_count.values()) <= max_conflict: # we are less than or equal the minium accepted number of conflict reads
        # calculate PS and HP
            for i in myvalues:
                # get the hp based on the hoghest number
                hp = max(hp_count, key = hp_count.get) # either 1 or 2
                ps =  most_frequent([i[1]  for i in myvalues if i[0] == hp])
                ps_dict[ps] = int(hp) if int(hp) > 1 else -1
        else: # Number of reads conflicting are higher than user suggestion
            ps = most_frequent([i[1] for i in myvalues])
            ps_dict[ps] = 0 # ['ps', 0]
    else:
        # Data are not conflict calculate normally:
        ps_dict = categorize_ps(myvalues)
    return ps_dict


def categorize_ps_up(myvalues: 'list', min_conflict: 'int') -> {}:
    myvalues = [i for i in myvalues if i is not None] # remove None
    ps_dict = {}
    if min_conflict > 0:
        min_conflict += 1
    for i in myvalues:
        hp = int(i[0])
        ps = i[1]
        if ps in ps_dict:
            if hp == 1:
                if ps_dict[ps] < 0:
                    ps_dict[ps] = ps_dict[ps] - 1
                else: #conflict
                    if min_conflict == 0:
                        ps_dict[ps] = 0
                    else:
                        min_conflict -= 1
                        ps_dict[ps] -= 1
            else: # means that it is haplotype 2 hp=2
                if ps_dict[ps] > 0:
                    ps_dict[ps] = ps_dict[ps] + 1
                else: #conflict
                    if min_conflict == 0:
                        ps_dict[ps] = 0
                    else:
                        min_conflict -= 1
                        ps_dict[ps] += 1
        else:
            if hp == 1:
                ps_dict[ps] = -1
            else:
                ps_dict[ps] = 1
    return ps_dict

def main():
    args = get_args()



if __name__ == "__main__":
    # del_test_840 = [['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['2', '531175'], ['1', '531175'], ['2', '531175'], ['2', '531175'], ['2', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175'], ['1', '531175'], ['2', '531175'], ['1', '531175']]
    #
    # hp = [i[0] for i in del_test_840]
    # hp_count = {i: hp.count(i) for i in hp}
    #
    # print(f'Number of reads supports HP {hp_count}')
    #
    # print(f'CONF 19 DEL 840 {categorize_ps_up(del_test_840, 19)}')
    # exit(1)
    main()
