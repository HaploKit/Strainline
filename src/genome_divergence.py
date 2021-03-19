import sys
import os

'''This program is used to compute the divergence between two genomes from two fasta files'''


def fasta_len(fa):
    sum_len = 0
    with open(fa, 'r') as fr:
        for line in fr:
            if line.startswith('>'):
                continue
            else:
                sum_len += len(line.strip())
    return sum_len


def cal_genome_divergence(fa1, fa2):
    paf_out = os.popen("minimap2 -cx asm20 -t 16 {} {} 2>/dev/null".format(fa1, fa2)).read().strip().split('\n')
    divergence = 1.0
    if not paf_out[0]:
        return divergence,divergence
    matched_len = 0  # identical bases
    ovlp_len = 0  # including mismatches and gaps
    for i, line in enumerate(paf_out):
        a = line.split('\t')
        matched_len += int(a[9])
        ovlp_len += int(a[10])
    fa1_len = fasta_len(fa1)
    fa2_len = fasta_len(fa2)

    fa1_uniq_ovlplen = os.popen("minimap2 -ax asm20 -t 16 {} {} 2>/dev/null|samtools sort - |samtools depth -|wc -l".
                                format(fa1, fa2)).read().strip()
    fa2_uniq_ovlplen = os.popen("minimap2 -ax asm20 -t 16 {} {} 2>/dev/null|samtools sort - |samtools depth -|wc -l".
                                format(fa2, fa1)).read().strip()
    fa1_uniq_ovlplen = int(fa1_uniq_ovlplen)
    fa2_uniq_ovlplen = int(fa2_uniq_ovlplen)
    global_divergence = round(1 - matched_len / (ovlp_len + fa1_len - fa1_uniq_ovlplen + fa2_len - fa2_uniq_ovlplen), 4)

    ovlp_divergence = round(1 - matched_len / ovlp_len, 4)  # only consider the overlap regions

    return global_divergence, ovlp_divergence


if __name__ == '__main__':
    fa1, fa2 = sys.argv[1:]
    # fa1='../data/1.fa'
    # fa2='../data/2.fa'
    global_divergence, ovlp_divergence = cal_genome_divergence(fa1, fa2)
    print('\n'+'-'*50)
    print('The global genome divergence is  : {}%'.format(global_divergence * 100))
    print('The overlap regions divergence is: {}%'.format(ovlp_divergence * 100))
    print('-'*50)

