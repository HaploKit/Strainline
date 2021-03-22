#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import argparse

INDEX = {
    "savage" : 1,
    "virus-vg" : 2,
    "vg-flow": 8,
    "abayesqr" : 5,
    "pehaplo" : 7,
    "predicthaplo" : 3,
    "shorah" : 4
}

FORMAT = ["method", "contig_id", "truth_id", "aln_len", "edit_dist",
            "mismatches", "ins_len", "del_len"]

def main():
    parser = argparse.ArgumentParser(prog='assembly_evaluation.py', description='Compare assembled contigs to ground truth haplotypes.')
    parser.add_argument('contigs', nargs='*', type=str)
    parser.add_argument('-a', '--assignments', dest='assignments', type=str, required=True)
    parser.add_argument('-n', '--num_strains', dest='num_strains', type=int, required=True)
    args = parser.parse_args()

    ref_count = args.num_strains
    dist_bins = range(0, 6)
    data = pd.read_table(args.assignments, names=FORMAT)
    # print(data)
    sns.set(style="ticks", context="paper", palette="muted")
    fig, axs = plt.subplots(1, 3, sharey=True, figsize=(10, 3), tight_layout=True)
    axs[0].set_title("Precision")
    axs[1].set_title("Recall")
    axs[2].set_title("F-measure")
    # print(data["method"].unique())
    for method in args.contigs:
        print(method)
        method_file = method.split('/')[-1]
        method_name = method_file.rstrip(".fasta")
        color = sns.color_palette()[INDEX[method_name]]
        ncontigs = fasta_len(method)
        if ncontigs == 0:
            continue
        subdata = data.loc[data["method"] == method_file]
        # print(subdata)
        stats = [compute_stats(subdata, dist, ref_count, ncontigs) for dist in dist_bins]
        precision = [x[0] for x in stats]
        recall = [x[1] for x in stats]
        f_score = [x[2] for x in stats]
        axs[0].plot(dist_bins, precision, label=method_name, color=color)
        axs[1].plot(dist_bins, recall, label=method_name, color=color)
        axs[2].plot(dist_bins, f_score, label=method_name, color=color)
    for ax in axs:
        ax.set_xticks(dist_bins)
        ax.set_xlabel("max % edit distance")
    plt.legend(loc='upper center', bbox_to_anchor=(1.35, 0.75), shadow=True, ncol=1)
    plt.savefig('prec_recall_fscore.eps')
    plt.show()
    return


def compute_stats(data, max_edit_perc, ref_count, ncontigs):
    true_positives = set()
    matched_ref = set()
    for index, record in data.iterrows():
        # print(record)
        edit_perc = record['edit_dist'] / record['aln_len'] * 100
        if edit_perc <= max_edit_perc:
            true_positives.add(record['contig_id'])
            matched_ref.add(record['truth_id'])
    recall = len(matched_ref) / ref_count
    precision = len(true_positives) / ncontigs if ncontigs > 0 else 0
    if precision + recall > 0:
        f_score = 2 * precision * recall /  (precision + recall)
    else:
        f_score = 0
    return precision, recall, f_score


def fasta_len(fname):
    with open(fname) as f:
        i = 0
        for line in f:
            if line[0] == '>':
                i += 1
    return i


if __name__ == '__main__':
    sys.exit(main())
