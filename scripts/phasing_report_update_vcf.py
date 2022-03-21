#!/usr/bin/env python3

"""
bamtools merge should be used before this script. where the vcf file should be merged with both paternal and maternal, respectively.
"""
import argparse
import sys, re

def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Produce phasing report',
                                     add_help=True, )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
    # parser.add_argument('input', help='Input file', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    # parser.add_argument('output', help='Output file', nargs="?", type=argparse.FileType('w'), default=sys.stdout)

    # parser.add_argument('input', nargs='?', help="Phased vcf file",
    #                          type=argparse.FileType('r'),
    #                          default=sys.stdin)
    # parser.add_argument('output', nargs='?', help="Output file if no file result will be directed to stander output",
    #                              type=argparse.FileType('w+'),
    #                              default=sys.stdout)
    parser.add_argument('-i', '--input', nargs='?', help="Phased vcf file", required=True)
    parser.add_argument('-o', '--output', nargs='?', help="Output file for blocks", required=True)
    parser.add_argument('-s', '--stat', nargs='?', help="Output statistics file for phased datat", required=True)
    parser.add_argument("-u", '--update_snps',  help="Output updated snp file", required=True)
    parser.add_argument('-t', '--tolerance', help="Percent of tolerance.", type=int, action='store', default=5)
    parser.add_argument('-n', '--min_snps', help="Minimum Number of SNPs per block.", type=int, action='store', default=10)

    parser.set_defaults(func=phase_filtering)
    args = parser.parse_args()
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()



def phase_filtering(args):
    # check if the input from stdin
    # if sys.stdin.isatty():
    #     if args.input.name.endswith("gz"):
    #         import gzip
    #         myfile = gzip.open(args.input.name, 'rt') # t is not a must normally it is default.
    #     else:
    #         myfile = args.input
    # else:
    myfile = args.input
    if myfile.endswith("gz"):
        import gzip
        myfile_open = gzip.open(myfile, 'rt', encoding='utf-8')
        #myfile_open = gzip.open(myfile, 'rt', encoding='utf-16')
    else:
        myfile_open = open(myfile, 'r')

    phasing_dictionary = {}


    with myfile_open as data_in, open(args.output, 'w') as data_out:
        nonphased_hetero = 0
        snp_number = 0
        maternal = 0
        paternal = 0
        homo_number = 0
        unknown_phased = 0
        non_sample_snp = 0
        for line in data_in:
            if line.startswith('#'):
                pass # chnage to print in output file
            else:
                snp_number += 1
                line_split = line.split()
                gt_flag, sample, father, mother = line_split[8:12]
                # gt_flag, sample, father, mother = line_split[8:12]
                if "1/1" in sample:
                    homo_number +=1
                # 9 is first sample followed by father and mother
                # ID f0/1 f1/0 m0/1 m1/0
                # paternal
                if "|" in sample and bool(re.search(r'\d',father)) and not bool(re.search(r'\d', mother)):
                    paternal += 1
                    gt_index = gt_flag.split(":").index("PS")
                    id = sample.split(":")[gt_index]

                    f0_1 = 0
                    f1_0 = 0
                    if sample[0] == "0":
                        f0_1 = 1
                    else:
                        f1_0 = 1

                    if id not in phasing_dictionary:
                        phasing_dictionary[id] = [f0_1, f1_0 , 0, 0]
                    else:
                        phasing_dictionary[id][0] += f0_1
                        phasing_dictionary[id][1] += f1_0
                # maternal
                elif "|" in sample and bool(re.search(r'\d', mother)) and not bool(re.search(r'\d',father)):
                    maternal += 1
                    gt_index = gt_flag.split(":").index("PS")
                    id = sample.split(":")[gt_index]

                    m0_1 = 0
                    m1_0 = 0
                    if sample[0] == "0":
                        m0_1 = 1
                    else:
                        m1_0 = 1

                    if id not in phasing_dictionary:
                        phasing_dictionary[id] = [0, 0, m0_1, m1_0]
                    else:
                        phasing_dictionary[id][2] += m0_1
                        phasing_dictionary[id][3] += m1_0
                # Unknown if it is right or wrong cause no equivliant in mother or father
                elif "|" in sample:
                    unknown_phased += 1
                elif "1/0" in sample or "0/1" in sample:
                #elif "1/1" not in sample and "." not in sample:
                    nonphased_hetero += 1 # is it hetero or homo zygot
                elif sample.startswith("."):
                    non_sample_snp += 1


        for k, v in phasing_dictionary.items():
            data_out.write("{}\t{}\n".format(str(k), "\t".join(map(str, v))))
        # print("Number SNPs: {snp}\nUnknown phased case: {unknown}\n \
        # Number of non-phased Hetero: {not_phased_hetero}\n \
        # Maternal phased: {mother}\nPaternal phased:  \
        # {father}\nDone".format(unknown=unknown_phased, not_phased_hetero=nonphased_hetero, snp=snp_number, mother=maternal, father=paternal))
        with open(args.stat, 'w') as stat_out:
            stat_out.write("\
            Number SNPs: {snp}\n\
            Homozygot number 1/1: {homo}\n\
            Unknown phased case: {unknown_cases}\n\
            Number of non-phased Hetero: {not_phased_hetero}\n\
            Total number of Phased SNPs: {total}\n\
            Maternal phased: {mother}\n\
            Paternal phased:{father}\n\
            SNP only in paternal: {no_snp}".format(unknown_cases=unknown_phased, homo=homo_number, not_phased_hetero=nonphased_hetero, snp=snp_number, mother=maternal, father=paternal, total=str(maternal+paternal), no_snp = non_sample_snp))


    # updating vcf phased snps
    chr = ""
    new_block_value = ""
    if args.update_snps:
        if myfile.endswith("gz"):
            import gzip
            myfile_open = gzip.open(myfile, 'rt', encoding='utf-8')
            #myfile_open = gzip.open(myfile, 'rt', encoding='utf-16')
        else:
            myfile_open = open(myfile, 'r')

        with myfile_open as data_in, open(args.update_snps , "w") as output:
            for line in data_in:
                if line.startswith('##'):
                    output.write(line)
                elif line.startswith("#"):
                    output.write("##INFO=<ID=parental-snps,Number=1,Type=Integer,Description=\"Total number of parental snps supporting the phsing block in called genotypes\">\n")
                    output.write("{}\n".format("\t".join(line.split()[:-2])))
                else:
                    line_split = line.split()
                    chr_value = line_split[0]
                    # identify wich chromosome we are using.
                    if chr_value != chr:
                        chr = chr_value
                        first_block = True
                    snp_format, format_value = line_split[8:10]
                    format_value_split = format_value.split(':')
                    if "PS" in snp_format and "|" in format_value.split(":")[0]:  # It is phased
                        # pritn(line)
                        # gt_value = format_value.split(":")[0]
                        block_value = format_value.split(":")[snp_format.split(":").index("PS")]
                        block_not_conflict = False
                        if block_value in phasing_dictionary:  # this snp have a similr one in parents vcf file.
                            block_not_conflict, gt_value = not_conflecting(phasing_dictionary[block_value], args)
                            # update the PS value.
                            if block_not_conflict:
                                # add +N of snps supoorting the block to p-snp
                                # udate the PS tag for each chromsome to be the first value in the first block (ps)
                                # write the updated line
                                # Update gt
                                if first_block:
                                    first_block = False
                                    new_block_value = block_value

                                format_value_split[0] = gt_value
                                # Update PS
                                format_value_split[snp_format.split(":").index("PS")] = new_block_value
                                line_split[9] = ":".join(format_value_split)
                                line_split[7] = line_split[7] + ";parental-snps=" + str(
                                    sum(phasing_dictionary[block_value]))
                                output.write("{}\n".format("\t".join(line_split[:-2])))
                                # print("block --> "+ new_block_value)
                            else:
                                # add -1 of snps supoorting the block to p-snp
                                # udate the PS tag for each chromsome to be the first value in the first block (ps)
                                # write the updated line
                                # Update gt
                                format_value_split[0] = gt_value
                                # Update PS
                                # format_value_split[snp_format.split(":").index("PS")] = new_block_value
                                line_split[9] = ":".join(format_value_split)
                                line_split[7] = line_split[7] + ";parental-snps=-" + str(
                                    sum(phasing_dictionary[block_value]))  # add sum
                                output.write("{}\n".format("\t".join(line_split[:-2])))
                                # print("block --> "+ new_block_value)
                        else:
                            # No information form parental about it
                            # keep it the same add 0 to p-snp flag
                            line_split[7] = line_split[7] + ";parental-snps=0"
                            output.write("{}\n".format("\t".join(line_split[:-2])))
                    elif bool(re.search(r'\d',format_value)):
                        # exit(format_value.split(":"))
                    # elif "/" in format_value.split(":")[0]:
                        # Write the line without change
                        line_split[7] = line_split[7] + ";parental-snps=."
                        output.write("{}\n".format("\t".join(line_split[:-2])))
                    # else:
                    #     print(line)





def hasNumbers(inputString):
     return any(char.isdigit() for char in inputString)

def is_not_conflict(block_snp_list, args):
    all_snps_in_block = sum(block_snp_list)
    if all_snps_in_block >= args.min_snps:
        # assuming that the list is formed like that F0|1 F1|0 M0|1 M1|0
        tolerance_percentage = args.tolerance * all_snps_in_block/100  # tolerance is 5%
        index_m, value_m = max(enumerate(block_snp_list[:2]), key=operator.itemgetter(1))
        index_f, value_f = max(enumerate(block_snp_list[2:]), key=operator.itemgetter(1))

        if  block_snp_list.count(0) == 3 or ( ( (index_m == 0 and index_f == 1) or (index_m == 1 and index_f == 0) ) and ( any(i <= 5*sum(block_snp_list[:2])/100 for i in block_snp_list[:2]) and any(i <= 5*sum(block_snp_list[2:])/100 for i in block_snp_list[2:])) ):
            return (True, "{}|{}".format(index_m, index_f))  # it means 0|1 or 1|0

def not_conflecting(block, args):
    tolerance = args.tolerance / sum(block) * 100
    max_index = block.index(max(block)) # bigest value in snps
    # assuming that the list is formed like that F0|1 F1|0 M0|1 M1|0
    # if sum(block) >= args.min_snps:
    if (max_index == 0 or max_index==3):
        non_conflict = block[0] + block[3]
        conflict = block[1] + block[2]
        gt = "0|1"  # parental|maternal
    elif (max_index == 1 or max_index ==2):
        non_conflict = block[1] + block[2]
        conflict = block[0] + block[3]
        gt = "1|0"
    return (conflict / (conflict + non_conflict) * 100 <= tolerance, gt)
    # else:
    #     return(False, "")




def main():
    args = get_args()



if __name__ == "__main__":
    main()
