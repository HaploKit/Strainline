'''
This script is used to sort long reads by read length and their overlaps,
and unify the strand which all reads are from.
'''
import os
import sys
import networkx as nx

from filter_ovlps import filter_ovlp


def sort_reads_by_len(in_fa, outdir, top_k):
    id2seq = {}
    id = ''

    with open(in_fa) as fr:
        for line in fr:
            if line.startswith('>'):
                id = line.strip()
            else:
                seq = line.strip()
                id2seq[id] = [seq, len(seq)]

    out_fa = outdir + '/' + '.'.join(os.path.basename(in_fa).split('.')[:-1]) + '.top' + str(top_k) + '.fa'
    fw = open(out_fa, 'w')
    i = 0
    for id, val in sorted(id2seq.items(), key=lambda d: d[1][1], reverse=True):
        i += 1
        if i <= top_k:
            fw.write(id + '\n' + val[0] + '\n')
    return out_fa


def cal_overlap(fasta, outdir, platform, threads, filter=True, min_ovlp_len=500, min_identity=0.0, o=200, r=0.8):
    '''
    calculate the overlaps of long reads
    '''
    paf = outdir + '/' + '.'.join(os.path.basename(fasta).split('.')[:-1]) + '.paf'
    filtered_paf = outdir + '/' + '.'.join(os.path.basename(fasta).split('.')[:-1]) + '.filtered.paf'
    if platform == 'pb':  # TODO: -c is necessary or not?
        # minimap = "minimap2 -x ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25  " + \
        #           "-t %s %s %s |cut -f 1-12 |fpa drop -i -m >%s" % (threads, fasta, fasta, paf)

        minimap = "minimap2 -x ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25  -t {} \
            {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}" \
            .format(threads, fasta, fasta, min_ovlp_len, min_identity, paf)

    elif platform == 'ont':
        # minimap = "minimap2 -x ava-ont -k15 -Xw5 -m100 -g10000 -r2000 --max-chain-skip 25  " + \
        #           "-t %s %s %s |cut -f 1-12 |fpa drop -i -m >%s" % (threads, fasta, fasta, paf)
        minimap = "minimap2 -x ava-ont -k15 -Xw5 -m100 -g10000 -r2000 --max-chain-skip 25 -t {} \
            {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}" \
            .format(threads, fasta, fasta, min_ovlp_len, min_identity, paf)
    else:
        raise Exception('wrong platform for reads type: pb/ont')
    os.system(minimap)
    if filter:
        filter_ovlp(paf, filtered_paf, min_ovlp_len, min_identity, o, r)
        return filtered_paf
    return paf


def count_lost_reads(fa, paf):
    num_reads = int(os.popen('cat {}|grep ">"|wc -l'.format(fa)).read())
    reads_in_ovlps = {}
    with open(paf) as fr:
        for line in fr:
            a = line.split()
            reads_in_ovlps[a[0]] = 1
            reads_in_ovlps[a[5]] = 1
    print('\n\nraw reads: {}, {} reads in filtered overlaps\n'.format(num_reads, len(reads_in_ovlps)))
    return num_reads, len(reads_in_ovlps)


def reverse_comp(seq):
    base2comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([base2comp[base] for base in rev_seq])
    return rev_comp_seq


def paf2graph(paf):
    edges = []
    with open(paf) as fr:
        for line in fr:
            a = line.strip().split()
            edges.append((a[0], a[5]))
    G = nx.Graph()
    G.add_edges_from(edges)
    return G


def sort_reads_by_ovlps(fa, paf, outdir, by_length):
    '''sort reads by length, keep the next read has overlap with the previous reads
    and unify the strand of all reads, assume '+'
    '''
    out_fa = outdir + '/' + '.'.join(os.path.basename(fa).split('.')[:-1]) + '.sorted.fa'
    if os.path.getsize(paf) == 0:
        os.system('cp {} {}'.format(fa,out_fa))
        return out_fa

    fw = open(out_fa, 'w')
    rr2strand = {}  # {read_i} {read_j} = '+'
    with open(paf, 'r') as fr:
        for line in fr:
            a = line.split()
            rr2strand.setdefault(a[0], {})[a[5]] = a[4]
            rr2strand.setdefault(a[5], {})[a[0]] = a[4]

    i = 0
    read2strand = {}
    used_reads = {}
    read2seq = {}
    first_read = ''
    with open(fa) as fr:
        read = ''
        for line in fr:
            if line.startswith('>'):
                i += 1
                read = line.strip().split()[0][1:]
                if i == 1:
                    first_read = read
                    read2strand[first_read] = '+'
            else:
                seq = line.strip()
                read2seq[read] = [seq, len(seq)]

    # output the first read info
    fw.write('>' + first_read + '\n' + read2seq[first_read][0] + '\n')

    used_reads[first_read] = 1
    reverse = 0  # count number of reads which need to reverse_complement
    read_i = 1
    next_read = first_read
    neighbors = {}
    neighbors[next_read] = 1
    while True:
        read_i += 1
        if read_i % 1000 == 0:
            print('processing the {} read...'.format(read_i))

        for r in rr2strand[next_read].keys():
            if r not in used_reads:
                neighbors[r] = 1
        del neighbors[next_read]
        # sort neighbor reads by length
        if by_length:
            if len(neighbors) > 0:
                next_read = sorted(neighbors.keys(), key=lambda r: read2seq[r][1], reverse=True)[0]
            else:
                raise Exception('No overlap found for read because the overlap graph is not connected !')
        else:
            for read in neighbors.keys():
                next_read = read
                break

        # check strand of reads
        for used_read in used_reads.keys():
            if used_read in rr2strand[next_read]:
                ovlp_strand = rr2strand[next_read][used_read]
                if read2strand[used_read] == '+':
                    if ovlp_strand == '+':
                        read2strand[next_read] = '+'
                    else:
                        read2strand[next_read] = '-'
                elif read2strand[used_read] == '-':
                    if ovlp_strand == '+':
                        read2strand[next_read] = '-'
                    else:
                        read2strand[next_read] = '+'
                break
                # TODO: check strand conflicts ?

        if read2strand[next_read] == '+':
            fw.write('>' + next_read + '\n' + read2seq[next_read][0] + '\n')
        elif read2strand[next_read] == '-':
            reverse += 1
            fw.write('>' + next_read + '\n' + reverse_comp(read2seq[next_read][0]) + '\n')
        else:
            raise Exception('Wrong read strand found: {}'.format(read2strand[next_read]))
        used_reads[next_read] = 1
        if len(used_reads) == len(rr2strand):  # all reads in overlap file are visited
            break
    fw.close()
    print('number of reverse_complement reads: {}'.format(reverse))
    return out_fa


if __name__ == '__main__':
    # fa = 'data/hiv/reads.fa'
    # outdir = 'data/hiv'
    # top_k=100
    fa, outdir, top_k, by_length, platform = sys.argv[1:]
    top_fa = sort_reads_by_len(fa, outdir, int(top_k))

    #mean read length=2.4kb, min_ovlp_len=400 in the first version
    paf = cal_overlap(top_fa, outdir, platform, threads=40, filter=True, min_ovlp_len=1000,
                      min_identity=0.1, o=400, r=0.8)
    count_lost_reads(top_fa, paf)

    # check if overlap graph is connected TODO: add a reference as the first read ?
    # G = paf2graph(paf)
    # if nx.is_connected(G):
    #     print('The overlap graph is connected.')
    # else:
    #     n = nx.algorithms.components.number_connected_components(G)
    #     print('The overlap graph is not connected and has {} components,\n'
    #           'Refiltering overlaps with more relaxed threshold is recommended.\n\n'.format(n))

    # run anyway
    sort_reads_by_ovlps(top_fa, paf, outdir, int(by_length))
