import sys


def reformat_fasta(in_fa, out_fa):
    id2seq = {}
    id = ''
    seq = ''
    with open(in_fa) as fr:
        for line in fr:
            if line.startswith('>'):
                if id:
                    id2seq[id] = [seq, len(seq)]
                    seq = ''
                id = line.strip()
            else:
                seq += line.strip()
    id2seq[id] = [seq, len(seq)]

    fw = open(out_fa, 'w')
    i = 0
    for id, val in sorted(id2seq.items(), key=lambda d: d[1][1], reverse=True):
        i += 1
        fw.write('>r' + str(i) + '\n' + val[0] + '\n')
    return out_fa


if __name__ == '__main__':
    in_fa, out_fa = sys.argv[1:]
    reformat_fasta(in_fa, out_fa)
