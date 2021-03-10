import sys
import os
from itertools import combinations
from multiprocessing import Pool

'''This program is used to remove redundant genomes based on sequence divergence'''


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
    ## Full genome/assembly alignment, intra-species asm-to-asm alignment
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
    fa_list_file,min_divergence,outdir = sys.argv[1:]
    min_divergence = float(min_divergence)
    fa_list=[]
    with open(fa_list_file) as fr:
        for line in fr:
            fa_list.append(line.strip())

    final_fastas={fa:1 for fa in fa_list} #the final fasta files after removing redundant genomes
    print(final_fastas)
    div_out = []
    p=0
    for fa1, fa2 in combinations(fa_list, 2):
        p+=1
        if p%100==0:
            print('the {} strain pair processing...'.format(p))
        global_divergence, local_divergence = cal_genome_divergence(fa1, fa2)
        div_out.append('\t'.join([fa1, fa2, str(global_divergence), str(local_divergence)]))
        if global_divergence < min_divergence:
            if fasta_len(fa1)>fasta_len(fa2):
                if fa2 in final_fastas:
                    del final_fastas[fa2]
            else:
                if fa1 in final_fastas:
                    del final_fastas[fa1]
    i=0
    with open(outdir+'strains.fa','w') as fw:
        for fa in final_fastas.keys():
            i+=1
            with open(fa) as fr:
                for line in fr:
                    if line.startswith('>'):
                        fw.write('>strain_{}'.format(i))
                    else:
                        fw.write(line)

    with open(outdir+'/genomes.divergence.txt', 'w') as fw:
        fw.write('\n'.join(div_out) + '\n')

