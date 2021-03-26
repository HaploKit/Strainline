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


def cal_genome_divergence(param):
    ## Full genome/assembly alignment, intra-species asm-to-asm alignment
    fa1, fa2,max_local_divergence = param
    paf_out = os.popen("minimap2 -cx asm20 -t 1 {} {} 2>/dev/null".format(fa1, fa2)).read().strip().split('\n')
    divergence = 1.0
    contained = 0  # if it is contained contig or not

    if not paf_out[0]:
        return (fa1, fa2, divergence, divergence,contained)

    matched_len = 0  # identical bases
    ovlp_len = 0  # including mismatches and gaps
    for i, line in enumerate(paf_out):
        a = line.split('\t')
        matched_len += int(a[9])
        ovlp_len += int(a[10])
    fa1_len = fasta_len(fa1)
    fa2_len = fasta_len(fa2)

    # number of bases which are covered in fa1 or fa2, which can be also obtained from paf file
    # this may be != fa1_len - ovlp_len (because of indels in overlap)
    # and != fa1_len - matched_len (because of mismatch bases not involved in matched_len)
    fa1_uniq_ovlplen = os.popen("minimap2 -ax asm20 -t 1 {} {} 2>/dev/null|samtools sort - |samtools depth -|wc -l".
                                format(fa1, fa2)).read().strip()
    fa2_uniq_ovlplen = os.popen("minimap2 -ax asm20 -t 1 {} {} 2>/dev/null|samtools sort - |samtools depth -|wc -l".
                                format(fa2, fa1)).read().strip()
    fa1_uniq_ovlplen = int(fa1_uniq_ovlplen)
    fa2_uniq_ovlplen = int(fa2_uniq_ovlplen)
    # global_divergence = 1 - matched / (similar+different)
    global_divergence = round(1 - matched_len / (ovlp_len + fa1_len - fa1_uniq_ovlplen + fa2_len - fa2_uniq_ovlplen), 4)
    local_divergence = round(1 - matched_len / ovlp_len, 4)  # only consider the overlap regions

    # Discard contained contigs
    fa1_oh = fa1_len - ovlp_len  # general overhang length of fa1
    fa2_oh = fa2_len - ovlp_len
    min_oh = 5
    # max_local_divergence = 0.001 #pacbio clr
    # max_local_divergence = 0.01 #for test, maybe for ont TODO
    if fa1_oh <= min_oh and local_divergence < max_local_divergence:
        contained = 1  # fa1 is contained contig
    elif fa2_oh <= min_oh and local_divergence < max_local_divergence:
        contained = 2

    return (fa1, fa2, global_divergence, local_divergence, contained)


if __name__ == '__main__':
    fa_list_file, max_global_divergence, outdir, threads,max_local_divergence = sys.argv[1:]
    max_global_divergence = float(max_global_divergence)
    max_local_divergence = float(max_local_divergence)
    fa_list = []
    with open(fa_list_file) as fr:
        for line in fr:
            fa_list.append(line.strip())

    final_fastas = {fa: 1 for fa in fa_list}  # the final fasta files after removing redundant genomes
    # print(final_fastas)
    div_out = []
    params = [(fa1, fa2,max_local_divergence) for fa1, fa2 in combinations(fa_list, 2)]

    pool = Pool(int(threads))
    res = pool.map(cal_genome_divergence, params, chunksize=1)  # ordered
    pool.close()
    pool.join()

    for line in res:
        # print('line:{}'.format(line))
        fa1, fa2, global_divergence, local_divergence, contained = line
        div_out.append('\t'.join([fa1, fa2, str(global_divergence), str(local_divergence), str(contained)]))
        if contained == 1:
            if fa1 in final_fastas:
                del final_fastas[fa1]
        elif contained == 2:
            if fa2 in final_fastas:
                del final_fastas[fa2]
        elif global_divergence < max_global_divergence:
            if fasta_len(fa1) > fasta_len(fa2):
                if fa2 in final_fastas:
                    del final_fastas[fa2]
            else:
                if fa1 in final_fastas:
                    del final_fastas[fa1]

    with open(outdir + '/haplotypes_divergence.txt', 'w') as fw:
        fw.write('\n'.join(div_out) + '\n')

    i = 1
    with open(outdir + '/haplotypes.fa', 'w') as fw:
        for fa in final_fastas.keys():
            with open(fa) as fr:
                for line in fr:
                    if line.startswith('>'):
                        continue
                    else:
                        if line.strip():
                            fw.write('>strain_{}\n'.format(i))
                            fw.write(line)
                            i += 1
