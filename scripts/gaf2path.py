import sys


def gaf2path(gaf, path_file, min_idt=0.5, min_match_ratio=0.5):
    path_info = {}
    num_records = 0
    with open(gaf, 'r') as fr:
        for line in fr:
            num_records += 1
            a = line.strip().split()
            name = a[0]
            query_len = int(a[1])
            match_len = int(a[3]) - int(a[2])
            strand = a[4]  # Strand relative to the path: "+" or "-"
            assert strand == '+'  # may exist '-' strand? TODO: if exist
            path_str = a[5]  # >4571>4573>4575 OR <662<661<659
            identity = float(a[14].split(':')[-1])
            out_cigar = '*'
            if path_str.startswith('>'):  # forward
                out_path = path_str[1:].replace('>', '+,')  # 1178+,1179+,1180+
                out_path = out_path + '+'  # add the last symbol
            else:
                out_path = '+,'.join(path_str.split('<')[::-1][:-1]) # reverse seq
                # out_path = path_str[1:].replace('<', '-,')
                out_path = out_path + '+'

            # only keep the record having the longest matched length if there are multiple alignments
            if (match_len / query_len >= min_match_ratio) and (identity >= min_idt):
                if name in path_info:
                    if match_len > path_info[name][0]:
                        path_info[name] = [match_len, '\t'.join(['P', name, out_path, out_cigar])]
                else:
                    path_info[name] = [match_len, '\t'.join(['P', name, out_path, out_cigar])]
    print('number of all alignments:{}'.format(num_records))
    print('number of retained alignments:{}'.format(len(path_info)))

    with open(path_file, 'w') as fw:
        fw.write('\n'.join([record[1] for _, record in path_info.items()]) + '\n')
    return


if __name__ == '__main__':
    gaf, path_file = sys.argv[1:]
    gaf2path(gaf, path_file)
