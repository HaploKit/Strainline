import sys


def path2seq(gfa, path_file, outfile):
    node2seq = {}
    with open(gfa, 'r') as fr:
        for line in fr:
            if line.startswith('S'):
                a = line.strip().split()
                node2seq[a[1]] = a[2]
    seq_list = []
    with open(path_file, 'r') as fr:
        for i, line in enumerate(fr):
            seq = '>' + str(i + 1) + '\n'
            for node in line.strip().split('+'):
                if node:
                    seq += node2seq[node]
            seq_list.append(seq)
    with open(outfile, 'w') as fw:
        fw.write('\n'.join(seq_list))
    return outfile


if __name__ == '__main__':
    gfa, path_file, outfile = sys.argv[1:]
    path2seq(gfa, path_file, outfile)
