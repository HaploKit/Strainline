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
    "abayesqr" : 5,
    "pehaplo" : 7,
    "predicthaplo" : 3,
    "shorah" : 4
}

FORMAT = ["method", "contig_id", "truth_id", "aln_len", "edit_dist",
            "mismatches", "ins_len", "del_len"]

def main():
    parser = argparse.ArgumentParser(prog='assembly_evaluation.py', description='Compare assembled contigs to ground truth haplotypes.')
    # parser.add_argument('contigs', nargs='*', type=str)
    parser.add_argument('-a', '--assignments', dest='assignments', type=str, required=True)
    parser.add_argument('-n', '--num_strains', dest='num_strains', type=int, required=True)
    parser.add_argument('-c', '--cov_list', dest='cov_list', type=str, required=True)
    args = parser.parse_args()

    ref_count = args.num_strains
    dist_bins = range(0, 6)
    data = pd.read_table(args.assignments, names=FORMAT)
    # print(data)
    sns.set(style="ticks", context="paper", palette="muted")
    # print(data["method"].unique())
    for cov in args.cov_list.split(','):
        print("{}x".format(cov))
        fig, axs = plt.subplots(1, 3, sharey=True, figsize=(10, 3), tight_layout=True)
        axs[0].set_title("Precision")
        axs[1].set_title("Recall")
        axs[2].set_title("F-measure")
        for method_name in INDEX.keys():
            print(method_name)
            color = sns.color_palette()[INDEX[method_name]]
            prec_list = []
            recall_list = []
            f_score_list = []
            for sample in range(1, 11):
                filename = "{}/sample{}.{}x.fasta".format(method_name, sample, cov)
                try:
                    ncontigs = int(fasta_len(filename))
                except FileNotFoundError as e:
                    continue
                if ncontigs == 0:
                    continue
                subdata = data.loc[data["method"] == filename]
                # print(subdata)
                stats = [compute_stats(subdata, dist, ref_count, ncontigs) for dist in dist_bins]
                prec_list.append([x[0] for x in stats])
                recall_list.append([x[1] for x in stats])
                f_score_list.append([x[2] for x in stats])
            nsamples = len(prec_list)
            if nsamples == 0:
                continue
            precision = np.mean(np.array(prec_list), axis=0)
            prec_std = np.std(np.array(prec_list), axis=0)
            recall = np.mean(np.array(recall_list), axis=0)
            recall_std = np.std(np.array(recall_list), axis=0)
            f_score = np.mean(np.array(f_score_list), axis=0)
            f_score_std = np.std(np.array(f_score_list), axis=0)
            # now plot these averages with error bars
            axs[0].errorbar(dist_bins, precision, yerr=prec_std, capsize=3,
                label=method_name, color=color)
            axs[1].errorbar(dist_bins, recall, yerr=recall_std, capsize=3,
                label=method_name, color=color)
            axs[2].errorbar(dist_bins, f_score, yerr=f_score_std, capsize=3,
                label=method_name, color=color)
        for ax in axs:
            ax.set_xticks(dist_bins)
            ax.set_xlabel("max % edit distance")
        plt.legend(loc='upper center', bbox_to_anchor=(1.35, 0.75), shadow=True, ncol=1)
        plt.savefig('prec_recall_fscore.{}x.eps'.format(cov))
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
