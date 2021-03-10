import sys
import os
import itertools
from multiprocessing import Pool

def get_posi_by_depth(lst, min_cov, w=5):
    '''
    get the start and end positions of fragments if there are low coverage regions
    in the middle of sequence, the input should be better trimmed at both ends.
    '''
    min_cov = min_cov  # * 0.8  # becuase of unstable coverage
    min_len = 500  # TODO, min length for output sequence
    k = int(len(lst) / w) + 1
    mean_cov = 0
    start = -1
    end = -1
    flag = False
    posi_list = []
    for i in range(k):
        mean_cov = sum([int(posi_cov.strip().split()[1]) for posi_cov in lst[i * w:(i + 1) * w]]) / w
        #         print('mean_cov:{}'.format(mean_cov))
        if mean_cov < min_cov:
            if start != -1:
                if (end - start) >= min_len:
                    posi_list.append((start, end, end - start))
                start = -1
                end = -1
                flag = False
            continue
        elif mean_cov >= min_cov and not flag:
            start = int(lst[i * w].strip().split()[0])
            end = int(lst[min(len(lst), (i + 1) * w) - 1].strip().split()[0])
            flag = True
        elif flag:

            end = int(lst[min(len(lst), (i + 1) * w) - 1].strip().split()[0])
    #             print('end:{}'.format(end))
    if flag and ((end - start) >= min_len):
        posi_list.append((start, end, end - start))  # add the last one
    return posi_list


def cov_based_misasm_removal_xx(param):
    i, reads_fa,outdir, min_cov, threads, trim_ends  = param
    contig_fa='contig.{}.fa'.format(i)
    os.system("mkdir -p {}".format(outdir))
    out_fa = outdir + '/xx.{}.fa'.format(i)

    # compute coverage for each base, filter secondary and supplementary alignments at first
    bam = outdir + "/" + str(i) + ".bam"
    os.system("minimap2 -a" + " --secondary=no -t 1 "  +" "+ contig_fa + " " + reads_fa + \
                  "  2>/dev/null |samtools view -hS -F 2048 -|samtools sort -@ 1 - >" + bam)
    # print("minimap2 -a" + " --secondary=no -t  " + str(threads) + contig_fa + " " + reads_fa + \
    #           "  2>/dev/null |samtools view -hS -F 2048 -|samtools sort -@ 24 - >" + bam)

    # TODO need to consider the low coverage region in the middle of sequence
    # @@ check why many errors in the middle of some super reads:
    # some reads are assigned into a wrong haplotype group or reads in the middle are belong to
    # one haplotype and reads in both ends are belong to the other, which may cause this mistake, therefore,
    # one should split the super reads into different parts and only keep the short fragments
    # (OR only keep the longest fragment) that are satisfied with min coverage requirement.
    """
    positions = os.popen("samtools depth " + bam + "|awk '$3>=" + str(min_cov) + \
                         "'|cut -f 2|sed -n '1p;$p'").read().strip()
    if positions:
        [start, end] = positions.split("\n")  # 1 based for sam
    else:
        print("{}.{}: no sequence satisfies the requirement of min coverage".format(i, hap))
        open(out_fa, 'w').close()  # new an empty file
        return
    """

    if os.path.getsize(bam) == 0:
        print("{}: no read can be aligned to the super read, skipping...".format(i),
              file=open("{}/log".format(outdir), 'a'))
        open(out_fa, 'w').close()  # new an empty file
        return
    else:
        positions = os.popen("samtools depth " + bam + "|awk '{print NR, $0}'" + "|awk '$4>=" + str(min_cov) + \
                             "'|cut -f 1 -d ' '|sed -n '1p;$p'").read().strip()
        if positions:
            [a, b] = [int(x) for x in positions.split("\n")]  # line number
            posi_cov_list = os.popen("samtools depth " + bam + "|cut -f 2,3").read().strip().split('\n')
            posi_list = get_posi_by_depth(posi_cov_list[(a - 1):b], min_cov, w=5)  # 1 based for sam
            print('posi_list:{}'.format(posi_list))
        else:
            print("contig.{}.fa: no sequence satisfies the requirement of min coverage".format(i),
                  file=open("{}/log".format(outdir), 'a'))
            open(out_fa, 'w').close()  # new an empty file
            return

    with open(contig_fa, "r") as fr:
        seq = fr.readline()
        seq = fr.readline().strip()

    sub_k = 1
    if trim_ends:
        for start, end, _ in posi_list:
            consensus = seq[(int(start) - 1):int(end)]
            if len(posi_list) == 1:
                head = ">c_{}\n".format(i)
            else:
                head = ">c_{}_sub{}\n".format(i, sub_k)
            with open(out_fa, 'a') as fw:
                fw.write(head + consensus + '\n')
            sub_k += 1
    else:
        # trim ends with low coverage and also break at the misassembled positions in the middle regions
        for start, end, _ in posi_list:
            if len(posi_list) == 1:
                consensus = seq[(int(start) - 1):int(end)]
                head = ">contig_{}\n".format(i)
            else:
                if sub_k == 1:
                    consensus = seq[:int(end)]
                elif sub_k == len(posi_list):
                    consensus = seq[(int(start) - 1):]
                else:
                    consensus = seq[(int(start) - 1):int(end)]
                head = ">contig_{}_sub{}\n".format(i, sub_k)
            with open(outdir+'/contig.{}.sub{}.fa'.format(i,sub_k), 'w') as fw:
                fw.write(head + consensus + '\n')
            sub_k += 1

    return


def cov_based_misasm_removal(param):
    i, reads_fa, all_bam, outdir, min_cov, threads, trim_ends = param
    # os.system("mkdir -p {}".format(outdir))
    out_fa = outdir + '/xx.{}.fa'.format(i)


    # compute coverage for each base, filter secondary and supplementary alignments at first
    bam = outdir + "/contig" + str(i) + ".bam"
    os.system("samtools view -b {} contig{} >{}".format(all_bam,str(i),bam))
    os.system("samtools index {}".format(bam))

    if os.path.getsize(bam) == 0:
        print("contig{}: bam is empty, skipping...".format(i),
              file=open("{}/log".format(outdir), 'a'))
        open(out_fa, 'w').close()  # new an empty file
        return
    else:
        positions = os.popen("samtools depth " + bam + "|awk '{print NR, $0}'" + "|awk '$4>=" + str(min_cov) + \
                             "'|cut -f 1 -d ' '|sed -n '1p;$p'").read().strip()
        if positions:
            [a, b] = [int(x) for x in positions.split("\n")]  # line number
            posi_cov_list = os.popen("samtools depth " + bam + "|cut -f 2,3").read().strip().split('\n')
            posi_list = get_posi_by_depth(posi_cov_list[(a - 1):b], min_cov, w=5)  # 1 based for sam
            print('# contig:{}, posi_list:{}'.format(i,posi_list))
        else:
            print("contig.{}.fa: no sequence satisfies the requirement of min coverage".format(i),
                  file=open("{}/log".format(outdir), 'a'))
            open(out_fa, 'w').close()  # new an empty file
            return

    contig="contig{}".format(i)
    seq = contig2seq[contig]

    sub_k = 1
    if trim_ends:
        for start, end, _ in posi_list:
            consensus = seq[(int(start) - 1):int(end)]
            if len(posi_list) == 1:
                head = ">c_{}\n".format(i)
            else:
                head = ">c_{}_sub{}\n".format(i, sub_k)
            with open(out_fa, 'a') as fw:
                fw.write(head + consensus + '\n')
            sub_k += 1
    else:
        # trim ends with low coverage and also break at the misassembled positions in the middle regions
        for start, end, _ in posi_list:
            if len(posi_list) == 1:
                consensus = seq[(int(start) - 1):int(end)]
                head = ">contig_{}\n".format(i)
            else:
                if sub_k == 1:
                    consensus = seq[:int(end)]
                elif sub_k == len(posi_list):
                    consensus = seq[(int(start) - 1):]
                else:
                    consensus = seq[(int(start) - 1):int(end)]
                head = ">contig_{}_sub{}\n".format(i, sub_k)
            with open(outdir + '/contig.{}.sub{}.fa'.format(i, sub_k), 'w') as fw:
                fw.write(head + consensus + '\n')
            sub_k += 1

    return


def clip_based_misasm_removal(read_fa, hap_fa, threads, outdir, min_clip_count):
    '''remove misassemblies based on clipped read alignments'''
    os.system("mkdir -p {}".format(outdir))
    sam = outdir + '/tmp.sam'
    min_clip_len = 5

    os.system(
        "minimap2 -a --secondary=no {} {}  -t {}|samtools view  -F 2048 - >{}".format(hap_fa, read_fa, threads, sam))
    hap2bpposi_tmp = {}
    hap2bpposi = {}

    with open(sam, 'r') as fr:
        for line in fr:
            a = line.split()
            hap, posi, cigar, seq = a[2], int(a[3]), a[5], a[9]
            read_len = len(seq)
            bp_posi = -1  # breakpoint position at target sequence
            # split cigar into list of numbers and characters all separately
            splitcigar = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]  # 25S100M2S
            if splitcigar[1] == 'S':  # soft-clipped, left
                if int(splitcigar[0]) < min_clip_len:
                    if splitcigar[-1] == 'S' and int(splitcigar[-2]) >= min_clip_len:
                        bp_posi = posi + read_len - int(splitcigar[-2])  # soft-clipped on both sides,such as 2S1000M50S
                    else:
                        continue
                else:
                    bp_posi = posi
            elif splitcigar[-1] == 'S':  # soft-clipped, right
                if int(splitcigar[-2]) < min_clip_len:
                    continue
                else:
                    bp_posi = posi + read_len - int(splitcigar[-2])
            else:
                continue
            if hap in hap2bpposi_tmp:
                hap2bpposi_tmp[hap].append(bp_posi)
            else:
                hap2bpposi_tmp[hap] = [bp_posi]

    for hap, posi_list in hap2bpposi_tmp.items():
        posi2count = {}
        for posi in posi_list:
            if posi in posi2count:
                posi2count[posi] += 1
            else:
                posi2count[posi] = 1
        for posi, count in posi2count.items():
            if count >= min_clip_count:
                print("## posi, count:{},{}".format(posi,count))
                if hap in hap2bpposi:
                    hap2bpposi[hap].append(posi)
                else:
                    hap2bpposi[hap] = [posi]

    hap2bpposi = {hap: sorted(posi_list) for hap, posi_list in hap2bpposi.items()}
    print('## haplotype to breakpoint position: {}'.format(hap2bpposi))
    hap=''
    out_info=[]
    i=0
    with open(hap_fa) as fr:
        for line in fr:
            sub_seqs=[]
            if line.startswith('>'):
                hap=line.strip()[1:]
            else:
                seq=line.strip()
                start=0
                if hap not in hap2bpposi:
                    i+=1
                    out_info.append('>hap_'+str(i)+'\n'+seq+'\n')
                    continue
                for bp_posi in hap2bpposi[hap]:
                    sub_seq = seq[start:(bp_posi-1)]
                    sub_seqs.append(sub_seq)
                    start=bp_posi-1

                sub_seqs.append(seq[start:]) #last sub sequence
                for s in sub_seqs:
                    if len(s)<500: #filter short sequences
                        continue
                    i+=1
                    out_info.append('>hap_' + str(i)+'\n'+s+'\n')

    with open(outdir+'/haplotypes.final.fasta','w') as fw:
        fw.write(''.join(out_info))
    return




if __name__ == '__main__':
    read_fa, contig_fa, outdir, threads= sys.argv[1:]
    os.system("mkdir -p {}".format(outdir))
    clip_based_misasm_removal(read_fa,contig_fa , threads, outdir, min_clip_count=10)
