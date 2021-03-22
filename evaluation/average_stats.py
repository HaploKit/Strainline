#!/usr/bin/env python3
import csv
import sys, os
import argparse

# df = pd.read_csv("test_quast_report.tsv", sep='\t')
# n_files = len(df.columns)-1
# print(n_files)


int_stats = ["# contigs", "N50", "NGA50"]
f1_stats = ["Genome fraction (%)"]
f2_stats = []
f3_stats = ["Error rate (%)", "Sensitivity", "PPV", "Rel freq error (%)",
    "Abs freq error (%)"]
output_stats = ["# contigs", "Genome fraction (%)", "N50", "NGA50",
    "Error rate (%)", "Sensitivity", "PPV", "Rel freq error (%)",
    "Abs freq error (%)"]

ER_stats = [
    "# N's per 100 kbp",
    "# mismatches per 100 kbp",
    "# indels per 100 kbp"
]
suppl_stats = [
    "# N's per 100 kbp",
    "# mismatches per 100 kbp",
    "# indels per 100 kbp",
    "Unaligned length"
]

def main():
    parser = argparse.ArgumentParser(prog='average_stats.py', description='Compute average assembly statistics from multiple input files.')
    parser.add_argument('reports', nargs='*', type=str)
    parser.add_argument('--dist', action='store_true')
    parser.add_argument('--skip_freq', action='store_true')
    parser.add_argument('--suppl', action='store_true')
    args = parser.parse_args()

    if not args.reports:
        print("No input files given. Use --help for usage information.")
        sys.exit(1)
    quast_reports = args.reports
    # quast_reports = ['test_quast_report.tsv']

    if args.dist:
        stats = []
        with open(quast_reports[0], 'r') as f:
            for line in f:
                strain = line.rstrip('\n').split('\t')[0]
                stats.append(strain)
        stats = sorted(stats[1:])
    else:
        if args.suppl:
            stats = output_stats + suppl_stats
        else:
            stats = output_stats
        if args.skip_freq:
            stats.remove("Rel freq error (%)")
            stats.remove("Abs freq error (%)")

    latex_table = []
    print('\nfilename\t' + '\t'.join(stats))
    for file in quast_reports:
        results = compute_average_stats(file, stats+ER_stats)
        result_line = file
        latex_line = file
        for stat in stats:
            if stat == "Error rate (%)":
                result = sum([results[stat] for stat in ER_stats])/1000
            else:
                try:
                    result = results[stat]
                except KeyError:
                    result = '-'
            # format output with desired floating point precision
            if result == '-':
                result_line += '\t{}'.format(result)
                latex_line += ' & {}'.format(result)
            elif stat in int_stats + suppl_stats:
                result_line += '\t{}'.format(int(round(result)))
                latex_line += ' & {}'.format(int(round(result)))
            elif stat in f1_stats:
                result_line += '\t{:.1f}'.format(result)
                latex_line += ' & {:.1f}'.format(result)
            elif stat in f2_stats:
                result_line += '\t{:.2f}'.format(result)
                latex_line += ' & {:.2f}'.format(result)
            elif stat in f3_stats:
                result_line += '\t{:.3f}'.format(result)
                latex_line += ' & {:.3f}'.format(result)
            elif args.dist:
                result_line += '\t{:.1f}'.format(result)
                latex_line += ' & {:.1f}'.format(result)
        print(result_line)
        latex_table.append(latex_line + '\\\\')
    print()
    print('\n'.join(latex_table))
    print()

    return


def compute_average_stats(file, stats="all"):
    results = {}
    with open(file, 'r') as tsv:
        tsv = csv.reader(tsv, delimiter='\t')
        for line in tsv:
            # print(line)
            stat = line[0]
            values = [x for x in line[1:] if x!='-']
            if stats=="all" or stat in stats:
                try:
                    av = sum([float(x) for x in values])/len(values)
                    results[stat] = av
                except ZeroDivisionError as e:
                    # all values were '-' i.e. not available
                    av = '-'
                # print(stat, av)
    return results

if __name__ == '__main__':
    sys.exit(main())
