#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import gzip

__author__ = "Vincent Luo"

usage = """%prog [options]

Split a fastq/fasta file into sub files randomly.
"""


def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument("--input", dest="input", type=str, required=True, help="input fastq/fasta filename")
    parser.add_argument("--output", dest="output", type=str, required=True, help="prefix for output files")
    parser.add_argument("--mode", dest="mode", type=str, required=True, help="fastq or fasta")
    parser.add_argument("--num", dest="num", type=int, required=True, help="number of reads in each bin")
    parser.add_argument("--pair", dest="pair", type=bool, required=False, help="if pair end reads, default: False")
    args = parser.parse_args()
    infile = args.input
    m = args.num  # number of reads in each bin
    mode = args.mode

    if infile.endswith('.gz'):
        n_reads = int(os.popen("zcat %s|wc -l" % (infile)).read().split()[0])
        # n_reads = int(os.popen("zless %s|wc -l" % (fastq)).read().split()[0]) #zcat does not work on Mac
    else:
        n_reads = int(os.popen("wc -l %s" % (infile)).read().split()[0])

    k = 0
    if mode == 'fastq':
        k = 8 if (args.pair) else 4
    elif mode == 'fasta':
        k = 4 if (args.pair) else 2
    else:
        print('ERROR: mode must be fastq or fasta')
        sys.exit(1)

    n_reads = int(n_reads / k)

    n_bins = int(n_reads / m)
    n_chosen_reads = n_bins * m

    if infile.endswith('.gz'):
        fq = gzip.open(infile, 'rb')
    else:
        fq = open(infile, 'r')

    i = 0
    item = ''
    id = 1
    id2item = {}

    for line in fq:
        i += 1
        if infile.endswith('.gz'):
            item += (line.decode().strip() + '\n')
        else:
            item += (line.strip() + '\n')
        if (i % k) == 0:
            seq_len = len(item.split('\n')[1])
            id2item[id] = [item, seq_len]
            item = ''
            id += 1
    fq.close()

    items = [t[0] for t in sorted(id2item.values(), key=lambda d: d[1], reverse=True)][:n_chosen_reads]

    random.seed(123)
    random.shuffle(items)

    for i in range(n_bins):
        with open("%s.sub_%s.fa" % (args.output, i + 1), 'w') as fw:
            fw.write(''.join(items[i * m:(i + 1) * m]))


if __name__ == '__main__':
    sys.exit(main())
