'''
Used to convert bam file into a hard/soft clipped fasta,
mainly for breaking chimeric reads or removing artifact sequences
in the both ends of reads.
The soft clipped bases will be removed in both primary and supplementary alignments.

Note: if the read is from - strand, the sequence in the bam is reversed and complementary,
but the order of CIGAR string is consistent with the sequence in the bam.

'''

import sys
import pysam
import itertools


def bam2fa(bam, fa, min_read_len):
    out = []
    with pysam.AlignmentFile(bam, 'rb') as fr:
        for line in fr:
            read_name, flag, _, _, _, cigar, _, _, _, seq = str(line).split()[:10]
            if flag=='4': #unmap reads
                continue
            # split cigar into list of numbers and characters all separately
            splitcigar = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
            # print(splitcigar)
            if splitcigar[1] == 'S':  # soft-clipped
                seq = seq[int(splitcigar[0]):]
            if splitcigar[-1] == 'S':
                seq = seq[:(-1 * int(splitcigar[-2]))]
            if len(seq) >= min_read_len:
                out.append('>' + read_name + '\n' + seq + '\n')

    with open(fa, 'w') as fw:
        fw.write(''.join(out))
    return

if __name__ == '__main__':
    # bam = '/Users/xiaoluo/Documents/CWI/project/vg/test.bam'
    # fa = 'xx.fa'
    # min_read_len = 300
    #
    bam, fa, min_read_len = sys.argv[1:]
    bam2fa(bam, fa, int(min_read_len))
