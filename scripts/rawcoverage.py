import argparse

import concurrent.futures as cf
import matplotlib.pyplot as plt
import numpy as np
# import multiprocessing as mp
import pandas as pd
import seaborn as sns

plt.switch_backend('agg')
from functools import partial
from os import path as opath
import ntpath
import sys
from Bio import SeqIO


def main():
    args = get_args()
    files = flat_list(args.input)
    nfiles = len(files)
    nworkers = min(nfiles, args.threads)
    # cpus  = mp.cpu_count()
    with cf.ProcessPoolExecutor(max_workers=nworkers) as executor, open(args.output, 'w') as data_out:
        # return pd.concat([i for i in executor.map(process_reads, files)], ignore_index=True)
        df = pd.concat([i for i in executor.map(process_reads, files)], ignore_index=True)
        data_out.write("Reads: {length}\n"
                       "Bases: {nbases}\n"
                       "Mean read length: {rmean}\n"
                       "Median: {rmdeian}\n"
                       "Max: {rmax}\n"
                       "Min: {rmin}\n"
                       "N50: {n50}". \
                       format(length=len(df), nbases=np.sum(df["lengths"]), rmean=np.mean(df["lengths"]),
                              rmdeian=np.median(df["lengths"]),
                              rmax=np.max(df["lengths"]),
                              rmin=np.min(df["lengths"]),
                              n50=get_N50(np.sort(df['lengths']))
                              ))

        plot_output = opath.join(opath.dirname(args.output), ntpath.basename(args.output).rsplit(".", 1)[0] + ".png")
        sns.set()
        myplot = sns.distplot(np.log(df['lengths']))
        myplot.set(xlabel='Log Read Length')
        myplot.get_figure().savefig(plot_output)

        # if i want to count reads c = s.groupby(['length']).size().reset_index(name='count')


def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for more info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Calulate statitics form fasta, fastq, fasta.gz and fastq.gz files ',
                                     add_help=True, )

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')

    parser.add_argument("-i", "--input", nargs="+",
                        help="<Required> one or more reads file ex: -i 1.fasta -i 2.fasta .... or -i 1.fasta  2.fasta",
                        action="append", required=True, metavar="FOO.fasta/q/gz")
    parser.add_argument("-o", "--output", help="<Required> output statitics file", metavar="FOO.txt")

    parser.add_argument("-t", "--threads", type=int, metavar='N', default=1,
                        help="<Optional> Number of threads default %(default)d")

    args = parser.parse_args()

    return args


def flat_list(my_list):
    """
    Transform list of lists to flat list
    :param my_list: list of lists ex: [[1],[1, 2], [a,v]]
    :return: [1, 1, 2, a, v]
    """
    return [element for each_list in my_list for element in each_list]


def process_reads(read_file):
    file_handle, file_type = open_handle(read_file)
    return (pd.DataFrame(
        data=[len(rec) for rec in SeqIO.parse(file_handle, file_type)],
        columns=["lengths"]).dropna())


def open_handle(myfile):
    if opath.isfile(myfile):
        if myfile.endswith('fastq.gz'):
            import gzip
            return gzip.open(myfile, 'rt'), "fastq"
        elif myfile.endswith('fasta.gz'):
            import gzip
            return gzip.open(myfile, 'rt'), "fasta"
        elif myfile.endswith('.fasta', ):
            return open(myfile, 'r'), 'fasta'
        elif myfile.endswith('.fastq'):
            return open(myfile, 'r'), 'fastq'
        # elif myfile.endswith("fastq.tar.gz"):
        #     import tarfile
        #     tar = tarfile.open(myfile, 'r:gz')#, 'fasta'
        #     for member in tar.getmembers():
        #          f = tar.extractfile(member)
        #          if f is not None:
        #              print(type(f))
        #              return open(f, 'r'), 'fastq'
        # elif myfile.endswith("fasta.tar.gz"):
        #     import tarfile
        #     tar = tarfile.open(myfile, 'r:gz')#, 'fasta'
        #     for member in tar.getmembers():
        #          f = tar.extractfile(member)
        #          if f is not None:
        #              return open(f, 'r'), 'fasta'
        else:
            sys.exit("This file {} is of unknow extesnion.".format(myfile))
    else:
        sys.exit("This file {} does not exist.".format(myfile))


def get_N50(read_lengths):
    return read_lengths[np.where(np.cumsum(read_lengths) >= 0.5 * np.sum(read_lengths))[0][0]]


if __name__ == "__main__":
    main()
