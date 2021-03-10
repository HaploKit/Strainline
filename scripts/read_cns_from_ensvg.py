#!/usr/bin/env python

## Generate consensus sequences of reads which are corrected by ensemble variation graphs ##

import os
import sys
from multiprocessing import Pool


def cns_from_ensvg(param):
    i, read, outdir = param
    if i % 1000 == 0:
        print('Processing the {} reads...'.format(i + 1))
    tmp_fa = outdir + '/' + read + '.fa'
    with open(tmp_fa, 'w') as fw:
        seqs = read2seq[read]
        for i in range(len(seqs)):
            fw.write(">{}\n{}".format(i, seqs[i]))

    consensus = os.popen('spoa -l 0 {}'.format(tmp_fa)).read()
    # consensus = os.popen('abpoa -m 1 {}'.format(cluster_fa)).read() #less computational cost

    consensus = consensus.strip().split('\n')[-1]
    read_item = ">{}\n{}\n".format(read, consensus)
    os.system('rm -rf {}'.format(tmp_fa))
    return read_item


if __name__ == '__main__':
    infile, outdir, threads = sys.argv[1:]
    read2seq = {}
    read = ''
    with open(infile) as fr:
        for i, line in enumerate(fr):
            if i % 2 == 0:
                read = line.strip()[1:]
            else:
                if read in read2seq:
                    read2seq[read].append(line)
                else:
                    read2seq[read] = [line]

    params = [(i, r, outdir) for i, r in enumerate(read2seq.keys())]

    pool = Pool(int(threads))
    res = pool.map(cns_from_ensvg, params, chunksize=1)  # ordered
    pool.close()
    pool.join()

    with open(outdir + '/reads.final.from_ensvg.fa', 'w') as fw:
        fw.write(''.join(res))
