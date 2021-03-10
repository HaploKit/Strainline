
'''convert MSA file to consensus sequence file'''

import sys
import re
from collections import defaultdict


msa_file = sys.argv[1]
col = int(sys.argv[2])  # which column represents the sequence
#
# msa_file='data/test.aln'
# col=1

pos_base2count = defaultdict(dict)
with open(msa_file) as fr:
    for line in fr:
        seq = line.strip().split()[col - 1]
        seq=re.sub('-+$','',seq) #trim '---' in the end
        flag=1
        for i, base in enumerate(seq):
            if base=='-' and flag: #skip '---' in the head
                continue
            flag=0
            if i in pos_base2count.keys():
                pos_base2count[i][base] += 1
            else:
                pos_base2count[i]['A'] = 0
                pos_base2count[i]['T'] = 0
                pos_base2count[i]['C'] = 0
                pos_base2count[i]['G'] = 0
                pos_base2count[i]['-'] = 0

consensus = ''
for pos in range(len(pos_base2count)):
    trust_base = ''
    flag = 0
    for base in pos_base2count[pos].keys():
        if pos_base2count[pos][base] >= flag:
            trust_base = base
            flag = pos_base2count[pos][base]
    if trust_base != '-':
        consensus += trust_base
print('>consensus\n' + consensus)
