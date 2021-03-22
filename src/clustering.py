import os
import sys
from multiprocessing import Pool
from filter_ovlps import filter_ovlp
from sort_reads import sort_reads_by_ovlps


def get_read2seq(file, mode='fastq'):
    read2seq = {}
    read = ''
    seq = ''
    if mode == 'fastq':
        with open(file) as fr:
            for i, line in enumerate(fr):
                if i % 4 == 0:
                    read = line.strip().strip('@')
                elif i % 4 == 1:
                    seq = line.strip()
                    read2seq[read] = seq
                else:
                    continue
    elif mode == 'fasta':
        with open(file) as fr:
            for i, line in enumerate(fr):
                if i % 2 == 0:
                    read = line.strip().strip('>')
                elif i % 2 == 1:
                    seq = line.strip()
                    read2seq[read] = seq
                else:
                    continue
    else:
        print("Error: unknown mode, only fastq or fasta permitted.")
        sys.exit(1)
    return read2seq


def write_fasta(k, reads, outdir, read2seq):
    fa = outdir + '/' + 'sread_cluster.' + str(k) + '.fa'
    with open(fa, 'w') as fw:
        for r in reads:
            fw.write(">" + r + "\n" + read2seq[r] + "\n")
    return fa


def sort_reads_by_len(in_fa, outdir):
    id2seq = {}
    id = ''

    with open(in_fa) as fr:
        for line in fr:
            if line.startswith('>'):
                id = line.strip()
            else:
                seq = line.strip()
                id2seq[id] = [seq, len(seq)]

    out_fa = outdir + '/' + '.'.join(os.path.basename(in_fa).split('.')[:-1]) + '.sort_by_len.fa'
    fw = open(out_fa, 'w')
    i = 0
    for id, val in sorted(id2seq.items(), key=lambda d: d[1][1], reverse=True):
        i += 1
        fw.write(id + '\n' + val[0] + '\n')
    return out_fa


def get_clusters(in_fa, outdir, topk, platform, threads, min_ovlp_len, min_identity, o, r, max_ovlps, min_sread_len,min_cluster_size):
    '''topk: top k seed read will be used
    '''
    id2fa = {}  # each read to its corresponding fasta
    id = 0
    info = ''
    with open(in_fa) as fr:
        for line in fr:
            if line.startswith('>'):
                info += line
            else:
                info += line
                id2fa[id] = info
                info = ''
                id += 1

    read2seq = get_read2seq(in_fa, 'fasta')
    k = 0
    used_reads = {}
    cluster_list = []
    for id in range(len(id2fa)):
        if k >= topk:
            break
        print('processing the {} seed read...'.format(k + 1))

        header, seq = id2fa[id].strip().split('\n')
        if len(seq) < min_sread_len:
            break
        sread = header[1:]  # seed read
        if sread in used_reads:
            continue
        with open(outdir + '/sread.fa', 'w') as fw:
            fw.write(header + '\n' + seq + '\n')
        # compute overlaps for seed read
        paf = outdir + '/sread.paf'
        filtered_paf = outdir + '/sread.filter.paf'

        # version 1
        os.system("minimap2 -cx ava-{}  -t {} {} {}|cut -f 1-12 >{}".format(platform, threads, outdir + '/sread.fa', in_fa, paf))

        # version 2:  do not use preset
        # os.system("minimap2 -c -k15 -w11  -t {} {} {}|cut -f 1-12 >{}".format(threads, outdir + '/sread.fa', in_fa, paf))
        # os.system("minimap2 -c -Hk23 -Xw7 -m100 -g10000 --max-chain-skip 25 -t {} {} {}|cut -f 1-12 >{}".format(threads, outdir + '/sread.fa', in_fa, paf))
        # os.system("minimap2 -cx asm5 -t {} {} {}|cut -f 1-12 >{}".format(threads, outdir + '/sread.fa', in_fa, paf)) # intra-species asm-to-asm alignment
        # os.system("minimap2 -X -c -k 21 -w 11 -s 60 -m 90 -n 2 -r 0 -A 4 -B 2 --end-bonus=100 -t {} {} {}|cut -f 1-12 >{}".format(threads, outdir + '/sread.fa', in_fa, paf)) # use parameters for short reads, from OGRE

        # if os.path.getsize(paf) <= 0:
        #     continue
        filter_ovlp(paf, filtered_paf, min_ovlp_len, min_identity, o, r, max_ovlps, rm_extra_ovlps=True)
        print('read overlaps filtering finished.')

        reads = []  # read IDs in this cluster

        if os.path.getsize(filtered_paf) == 0:
            reads.append(sread)
        else:
            with open(filtered_paf, 'r') as fr:
                for line in fr:
                    a = line.split()
                    reads.extend([a[0], a[5]])
        reads = set(reads)
        if len(reads) < min_cluster_size:  # iter2,only one contig is also fine
            # if len(reads)<5:
            continue
        sread_cluster_fa = write_fasta(k, reads, outdir, read2seq)
        for read in reads:
            used_reads[read] = 1

        # sort reads
        sorted_fa = sort_reads_by_ovlps(sread_cluster_fa, filtered_paf, outdir, by_length=1)
        cluster_list.append(sorted_fa)
        k += 1
    return cluster_list


def get_consensus(param):
    # compute consensus
    i, cluster_fa, outdir = param

    #sort reads before run spoa, already sorted in get_clusters()
    consensus = os.popen('spoa -l 0 {}'.format(cluster_fa)).read()
    # consensus = os.popen('abpoa -m 1 {}'.format(cluster_fa)).read() #less computational cost,but result seems not good as spoa

    consensus = consensus.strip().split('\n')[-1]
    header = os.popen('cat {}|head -1'.format(cluster_fa)).read()
    out_fa = outdir + '/contig.' + str(i) + '.fa'
    with open(out_fa, 'w') as fw:
        fw.write(header + consensus + '\n')
    return out_fa


def get_consensus_parallel(cluster_list, threads, outdir):
    pool = Pool(threads)
    params = [(i, cluster, outdir) for i, cluster in enumerate(cluster_list)]
    pool.map(get_consensus, params, chunksize=1)  # ordered
    pool.close()
    pool.join()
    # with open(outdir + '/contigs.fa', 'w'd) as fw:
    #     fw.write(''.join(out))
    return


if __name__ == '__main__':
    # topk: top k seed reads
    in_fa, outdir, topk, platform, threads, min_ovlp_len, min_identity, o, r, max_ovlps, min_sread_len = sys.argv[1:]
    cluster_list = get_clusters(in_fa, outdir, int(topk), platform, int(threads), int(min_ovlp_len),
                                float(min_identity), int(o), float(r), int(max_ovlps), int(min_sread_len),min_cluster_size=1)
    get_consensus_parallel(cluster_list, int(threads), outdir)
    os.system('rm -rf {}/sread*'.format(outdir))
