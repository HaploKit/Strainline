#!/usr/bin/env python3
import os
import sys
import subprocess
import copy
import multiprocessing as mp
import numpy as np # array, arange, minimum, add (levenshtein)
import itertools # groupby()
import argparse

def main():
    parser = argparse.ArgumentParser(prog='assembly_evaluation.py', description='Compare assembled contigs to ground truth haplotypes.')
    parser.add_argument('-c', '--contigs', type=str, required=True)
    parser.add_argument('-gt', '--ground_truth', dest='truth', type=str, required=True)
    parser.add_argument('-m', '--max_edit', dest='max_edit', type=float, default=0)
    parser.add_argument('-o', '--outdir', dest='outdir', default=".", help="write output files to this directory")
    parser.add_argument('-f', '--freq_truth', dest='freq_truth', help="file containing true frequency per ground truth strain")
    parser.add_argument('--total_cov', dest='total_cov', type=float, required=True)
    args = parser.parse_args()

    print("-----------------------------------------------")
    if not args.contigs:
        print("No contig file given. Aborting evaluation.")
        sys.exit(1)
    else:
        print("Contig file = {}".format(args.contigs))
        print("Truth file = {}".format(args.truth))

    # read contig abundance estimates and true haplotypes
    contigs, ab_est = read_fasta(args.contigs, read_ab=True)
    ground_truth, empty = read_fasta(args.truth)
    freq_truth = read_frequencies(args.freq_truth)
    print("Number of contigs: {}".format(len(contigs)))
    print("Number of reference sequences: {}".format(len(ground_truth)))
    if not contigs:
        print("No contigs found, exiting.")
        sys.exit(1)

    # assign contigs to strains -- note that:
    # (1) a contig may be assigned to multiple strains (conserved region)
    # (2) multiple contigs could map to the same strain, in the same region,
    # due to uncorrected errors. To avoid this situation we only evaluate
    # perfectly matching contigs
    # format: {contig id : [score, aln, aln_range]}
    contig_assignments = assign_contigs(args.contigs, args.truth, args.max_edit,
                                        args.outdir)

    # how many contigs match perfectly?
    evaluation_ratio = len(contig_assignments) / len(contigs)
    discarded = len(contigs) - len(contig_assignments)
    print("# contigs excluded from evaluation = {}".format(discarded))
    print("Evaluation ratio = {:.2f}\n".format(evaluation_ratio))

    total_assembly_freq = 0
    assembled_seqs = set([x[1][1] for l in contig_assignments.values() for x in l])
    for seq, freq in freq_truth.items():
        if seq in assembled_seqs:
            total_assembly_freq += freq

    # calculate true contig abundances
    results = []
    with open("{}/contig_abundance_estimates.txt".format(args.outdir), 'w') as f:
        for contig_id, assignment in contig_assignments.items():
            true_freq = calculate_freq(assignment, freq_truth, total_assembly_freq)
            true_ab = true_freq * args.total_cov
            estimate = ab_est[contig_id]
            results.append((true_ab, estimate))
            pos = assignment[0][1][2]
            f.write("{} {} {}\n".format(true_ab, estimate, pos))
            print(contig_id, assignment[0][1][1], len(contigs[contig_id]),
                true_ab, estimate)
    if not results:
        print("No results calculated, exiting.")
        sys.exit(1)

    # evaluate estimates
    total_error = 0
    abs_err = 0
    for (v1, v2) in results:
        if v1 + v2 > 0:
            total_error += abs(v1-v2) / (0.5 * (v1+v2))
            abs_err += abs(v1-v2) / args.total_cov
    rel_error = total_error / len(results) * 100
    abs_err = abs_err / len(results) * 100
    print("Relative abundance estimation error = {:.1f}%".format(rel_error))
    # print("Absolute abundance estimation error = {:.1f}%".format(abs_err))
    print("-----------------------------------------------")
    return


def calculate_freq(assignment, freq_truth, total_assembly_freq):
    true_freq = 0
    for aln_info in assignment:
        [score, aln, aln_range] = aln_info
        truth_id = aln[1]
        true_freq += freq_truth[truth_id]
    if true_freq == 0:
        print(assignment)
    assert true_freq > 0
    return true_freq / total_assembly_freq


def assign_contigs(contig_file, truth_file, max_edit, outdir):
    contigs, empty = read_fasta(contig_file)
    ground_truth, empty = read_fasta(truth_file)
    contigs2aln= {contig_id : [] for contig_id in contigs.keys()}

    # run bwa mem -a to align contigs to ground truth
    print("\nAligning contigs against ground truth...")
    sam_file = "{}/contigs_to_truth.sam".format(outdir.rstrip('/'))
    subprocess.check_call("bwa index {} 2>/dev/null".format(truth_file),
        shell=True)
    subprocess.check_call(
        "bwa mem -a -L1000 -t 12 {0} {1} > {2} 2>/dev/null".format(
            truth_file, contig_file, sam_file),
        shell=True)
    print()
    # score all alignments per contig and select the good alignment(s)
    [alignments, unmapped_ids] = read_sam(sam_file)
    print("#alignments: {}".format(len(alignments)))
    print("#unmapped: {}".format(len(unmapped_ids)))
    # score and store alignments
    for aln in alignments:
        #print aln
        [contig_id, ref_id] = aln[0:2]
        contig_seq = contigs[contig_id]
        truth_seq = ground_truth[ref_id]
        [score, aln_range] = score_alignment(aln, contig_seq, truth_seq)
        if score[0] <= max_edit:
            contigs2aln[contig_id].append([score, aln, aln_range])
    os.remove(sam_file)

    # select optimal alignment(s)
    contig_assignments = {}
    for contig, aln_list in contigs2aln.items():
        min_score = 1
        opt_aln = []
        #print contig
        #print aln_list
        for i in range(len(aln_list)):
            info1 = aln_list[i]
            if info1 in opt_aln:
                # alignment already paired with another alignment
                continue
            scores1 = info1[0]
            combined = False
            for j in range(i+1, len(aln_list)):
                info2 = aln_list[j]
                # check if these alignments can be combined
                if info2[1][1] != info1[1][1]: # not to same ref sequence
                    continue
                elif not (info2[2][1] <= info1[2][0]
                            or info2[2][0] >= info1[2][1]):
                    # overlapping alignments
                    if info2[2][1] >= info1[2][0] and info1[2][1] >= info2[2][0]:
                        ov_len = info1[2][1]-info2[2][0]
                    elif info1[2][1] >= info2[2][0] and info2[2][1] >= info1[2][0]:
                        ov_len = info2[2][1]-info1[2][0]
                    else:
                        print("no overlap???")
                        sys.exit(1)
                    print("\noverlap of length {} for contig {} on ref {}\n".format(ov_len, contig, info1[1][1]))
                    ref_len = info1[2][2]
                    if ov_len > 0.1*ref_len:
                        print("overlap too long, skipping")
                        continue
                # compute combined score (assigning contig to 2 allele sequences)
                combined = True
                combined_score = (info1[0][0]*info1[2][2]
                        + info2[0][0]*info2[2][2]) / (info1[2][2] + info2[2][2])
                # now check if this alignment is optimal
                if combined_score == min_score:
                    ref_id = info1[1][1]
                    new_aln = True
                    for aln in opt_aln:
                        if aln[1][1] != ref_id:
                            continue
                        else:
                            new_aln = False
                            break
                    if new_aln:
                        opt_aln.append(info1)
                        opt_aln.append(info2)
                elif combined_score < min_score:
                    opt_aln = [info1, info2]
                    min_score = combined_score
            # if there is no pairing, check single alignemnt
            if not combined:
                scores = scores1
                if scores[0] == min_score:
                    ref_id = info1[1][1]
                    new_aln = True
                    for aln in opt_aln:
                        if aln[1][1] != ref_id:
                            continue
                        else:
                            new_aln = False
                            break
                    if new_aln:
                        opt_aln.append(info1)
                elif scores[0] < min_score:
                    opt_aln = [info1]
                    min_score = scores[0]
        if opt_aln:
            contig_assignments[contig] = opt_aln
        count = len(opt_aln)
        if count > 1:
            print("NOTE: contig {} has {} assignments!".format(contig, count))
    return contig_assignments


def score_alignment(aln, contig, truth):
    [seq_id, ref_id, pos, ori, cigar] = aln
    if ori == "-":
        contig = revcomp(contig)
    # split cigar into numbers and characters all separately
    splitcigar = ["".join(x) for _, x in itertools.groupby(cigar,
                    key=str.isdigit)]
    # count mismatches, insertions, deletions and N's
    mismatch_count = 0
    ins_count = 0
    ins_len = 0
    del_count = 0
    del_len = 0
    N_count = 0
    length = len(contig)
    # keep track of position and cigar index
    contig_pos = 0
    truth_pos = pos
    start_pos = pos
    i = 0
    while i+1 < len(splitcigar):
        aln_len = int(splitcigar[i])
        aln_type = splitcigar[i+1]
        if aln_type == "S" or aln_type == "H":
            if i == 0: # front end clipped
                # if pos > 0:
                #     del_count += 1
                #     del_len += pos
                clipped_truth = min(aln_len, pos)
                length -= aln_len - clipped_truth # correct for overhanging contig length
                start_pos -= clipped_truth
                contig_pos += aln_len
            elif i > 0: # back end clipped
                clipped_truth = min(aln_len, len(truth)-truth_pos)
                # clipped_truth = len(truth) - truth_pos
                # if clipped_truth < 0:
                #     del_count += 1
                #     del_len += clipped_truth
                length -= aln_len - clipped_truth
                contig_pos += aln_len
                # truth_pos = len(truth) - clipped_truth
                truth_pos += clipped_truth
            # compare (partially) aligned sequences
            contig_part = contig[contig_pos-clipped_truth : contig_pos]
            truth_part = truth[truth_pos-clipped_truth : truth_pos]
            [mismatches, Ns] = count_mismatches(contig_part, truth_part)
            mismatch_count += mismatches
            N_count += Ns
        elif aln_type == "I":
            if i > 0 and i < len(splitcigar)-2:
                ins_count += 1
                ins_len += aln_len
            N_count += contig[contig_pos:contig_pos+aln_len].count('N')
            contig_pos += aln_len
        elif aln_type == "D":
            del_count += 1
            del_len += aln_len
            truth_pos += aln_len
        elif aln_type == "M":
            if contig_pos + aln_len > len(contig):
                print("contig {} too short?".format(seq_id))
                print(contig_pos, aln_len, len(contig))
                print(pos, cigar)
            if truth_pos + aln_len > len(truth):
                print("truth too short?")
                print(truth_pos, aln_len, len(truth))
                print(pos, cigar)
            contig_part = contig[contig_pos : contig_pos + aln_len]
            truth_part = truth[truth_pos : truth_pos + aln_len]
            [mismatches, Ns] = count_mismatches(contig_part, truth_part)
            mismatch_count += mismatches
            N_count += Ns
            truth_pos += aln_len
            contig_pos += aln_len
        else:
            print("ERROR: cigar string not recognized.")
            sys.exit(1)
        i += 2
    # compute percent identity
    if length > 0:
        score = (mismatch_count + ins_len + del_len) / float(length)
    else:
        score = 0
    aln_scores = [score, mismatch_count, ins_count, ins_len, del_count, del_len,
                    N_count]
    assert start_pos >= 0
    assert truth_pos <= len(truth)
    assert start_pos <= truth_pos
    aln_range = [start_pos, truth_pos, len(truth)]
    edit_distance = mismatch_count + ins_len + del_len
    return [aln_scores, aln_range]

def count_mismatches(seq1, seq2):
    assert len(seq1) == len(seq2)
    mismatches = 0
    Ns = 0
    for i in range(len(seq1)):
        if seq1[i] == "N": # truth has N
            Ns += 1
        elif seq2[i] == "N": # contig has N
            Ns += 1
        elif seq1[i] != seq2[i]:
            mismatches += 1
    return [mismatches, Ns]


def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    revcomp = "".join(complement.get(base, base) for base in reversed(seq))
    assert len(seq) == len(revcomp)
    return revcomp


def read_frequencies(freq_truth):
    # read true frequencies
    true_frequencies = {}
    with open(freq_truth, 'r') as f:
        for line in f:
            [seq_id, freq] = line.strip('\n').split()
            true_frequencies[seq_id] = float(freq)
    return true_frequencies


def read_fasta(filename, read_ab=False):
    # returns ID to sequence dict
    id2seq = {}
    ab_est = {}
    with open(filename, 'r') as f:
        seq_id = ""
        seq = ""
        for line in f:
            if line[0] == '>':
                if seq_id != "" and seq != "":
                    id2seq[seq_id] = seq
                    if read_ab:
                        ab_est[seq_id] = ab # estimated abundance
                seq_id = line.lstrip('>').rstrip('\n').split()[0]
                if read_ab:
                    try:
                        ab = float(line.lstrip('>').rstrip('\n').split()[-1].lstrip('ab=').rstrip('x'))
                    except ValueError:
                        print("WARNING: could not read abundance estimates from fasta")
                        read_ab = False
                seq = ""
            else:
                seq += line.rstrip('\n')
        # add final entry
        if seq_id != "" and seq != "":
            id2seq[seq_id] = seq
            if read_ab:
                ab_est[seq_id] = ab # estimated abundance
    return id2seq, ab_est


def read_sam(filename):
    # returns a list of all alignments, where each alignment is presented in the
    # following format: [seq_id, ref_id, pos, ori, cigar]
    aln_list = []
    unmapped_ids = []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == "@":
                continue
            line = line.rstrip('\n').split()
            [seq_id, flag, ref_id, pos] = line[0:4]
            cigar = line[5]
            bits_flag = power_find(int(flag))
            if 4 in bits_flag: # unmapped contig
                unmapped_ids.append(seq_id)
                continue
            elif 16 in bits_flag: # reversed contig
                ori = "-"
            else:
                ori = "+"
            aln_list.append([seq_id, ref_id, int(pos)-1, ori, cigar]) # SAM-format is 1-based; switching to 0-based here
    return [aln_list, unmapped_ids]

def power_find(n):
    try:
        int(n)
    except TypeError:
        print("power_find TypeError")
        print("n = {}".format(n))
        sys.exit(1)
    result = []
    binary = bin(n)[:1:-1]
    for x in range(len(binary)):
        if int(binary[x]):
            result.append(2**x)
    return result


if __name__ == '__main__':
    sys.exit(main())