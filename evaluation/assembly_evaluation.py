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
    parser.add_argument('contigs', nargs='*', type=str)
    parser.add_argument('-gt', '--ground_truth', dest='truth', type=str, required=True)
    parser.add_argument('-m', '--max_edit', dest='max_edit', type=float, default=0.05)
    parser.add_argument('-c', '--min_HC', dest='min_HC', type=float, default=0.5)
    parser.add_argument('-o', '--outdir', dest='outdir', default=".", help="write output files to this directory")
    parser.add_argument('-f', '--freq_truth', dest='freq_truth', help="file containing true frequency per ground truth strain")
    args = parser.parse_args()

    if not args.contigs:
        print("No contig files given. Aborting evaluation.")
        sys.exit(1)

    main_outfile = args.outdir.rstrip('/') + '/extra_stats.tsv'
    dist_outfile = args.outdir.rstrip('/') + '/min_dist.tsv'
    freq_outfile = args.outdir.rstrip('/') + '/freq_true_est.tsv'
    assignment_outfile = args.outdir.rstrip('/') + '/contig_assignments.tsv'

    if os.path.exists(freq_outfile):
        os.remove(freq_outfile)
    if os.path.exists(assignment_outfile):
        os.remove(assignment_outfile)

    truth_file = args.truth # "/home/jasmijn/HLA_real/ground_truth/NA12878.HLA-A.fasta"
    contig_files = args.contigs # "/home/jasmijn/HLA_real/GIAB_NA12878/savage-lc/contigs_diploid.fasta"
    bwa_mode = True

    # read in reference sequences
    ground_truth, empty = read_fasta(truth_file)

    # process assemblies
    header_line = "Assembly"
    sens_line = "Sensitivity"
    ppv_line = "PPV"
    freq_line = "Rel freq error (%)"
    freq_line2 = "Abs freq error (%)"
    min_dist_data = {}
    for file in contig_files:
        header_line += "\t{}".format(file)
        stats = evaluate_assembly(file, ground_truth, bwa_mode, args.max_edit,
                    args.freq_truth, truth_file, args.min_HC, freq_outfile,
                    args.outdir, assignment_outfile)
        sens_line += "\t{:.2f}".format(stats[0])
        ppv_line += "\t{:.2f}".format(stats[1])
        freq = stats[2]
        freq_abs = stats[3]
        if freq >= 0:
            freq_line += "\t{:.1f}".format(freq)
            freq_line2 += "\t{:.1f}".format(freq_abs)
        else:
            freq_line += "\t-"
            freq_line2 += "\t-"
        min_dist_data[file] = stats[4]

    with open(main_outfile, 'w') as f:
        f.write(header_line + '\n')
        f.write(sens_line + '\n')
        f.write(ppv_line + '\n')
        f.write(freq_line + '\n')
        f.write(freq_line2 + '\n')

    with open(dist_outfile, 'w') as f:
        f.write(header_line + '\n')
        for strain in min_dist_data[contig_files[0]].keys():
            line = strain
            for file in contig_files:
                dist = min_dist_data[file][strain]
                if dist == 100:
                    line += "\t-"
                else:
                    line += "\t{:.1f}".format(min_dist_data[file][strain])
            f.write(line + '\n')

    return

################################################################################

def evaluate_assembly(contig_file, ground_truth, bwa_mode, max_edit, freq_truth,
            truth_file, min_HC, freq_out, outdir, assignment_out):
    contigs, ab_est = read_fasta(contig_file, True)
    print("------------------------------------------")
    print("Number of contigs: {}".format(len(contigs)))

    contigs2aln= {}

    if bwa_mode:
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
        for contig in contigs:
            contigs2aln[contig] = []
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
    else:
        print("\nPerforming pairwise alignments...")
        for contig_id, contig_seq in contigs.items():
            aln_list = []
            for truth_id, truth_seq in ground_truth.items():
                for ori in ["+", "-"]:
                    if ori == "+":
                        dt, bt = align2(truth_seq, contig_seq)
                    else:
                        dt, bt = align2(truth_seq, revcomp(contig_seq))
                    cigar, pos = backtrace_to_cigar(bt)
                    aln = [contig_id, truth_id, pos, ori, cigar]
                    [score, aln_range] = score_alignment(aln, contig_seq, truth_seq)
                    if score[0] <= max_edit:
                        aln_list.append([score, aln, aln_range])
            contigs2aln[contig_id] = aln_list

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
        contig_assignments[contig] = opt_aln
        count = len(opt_aln)
        if count > 1:
            print("NOTE: contig {} has {} assignments!".format(contig, count))

    # build a dict mapping ground truth fragments to contig alignments
    truth2aln = {}
    for ref_id in ground_truth:
        truth2aln[ref_id] = []
    for contig, aln_list in contig_assignments.items():
        for info in aln_list:
            scores = info[0]
            aln = info[1]
            ref_id = aln[1]
            aln_range = info[2]
            truth2aln[ref_id].append([aln, scores, aln_range])

    # for ref, aln_list in truth2aln.items():
    #     print ref
    #     for aln in aln_list:
    #         print aln

    # evaluate assembly
    contig_lengths = [len(seq) for contig,seq in contigs.items()]
    total_assembly_len = sum(contig_lengths)
    N50_all_contigs = compute_NX(contig_lengths, 50)
    outfile_line, sensitivity, ppv, min_dist_map = get_assembly_stats(truth2aln,
        total_assembly_len, N50_all_contigs, ground_truth, len(contigs2aln),
        assignment_out, contig_file)
    if freq_truth and len(ab_est) > 0:
        av_rel_err, av_abs_err, median = check_frequencies(
            freq_truth, ab_est, truth2aln, contigs, ground_truth, min_HC,
            freq_out)
        if av_rel_err >= 0:
            print("average relative abundance error: {:.1f}%".format(av_rel_err))
            print("average absolute abundance error: {:.1f}%".format(av_abs_err))
            print("median abundance error: {:.1f}%\n".format(100*median))
    else:
        av_rel_err = -1
        av_abs_err = -1

    return sensitivity, ppv, av_rel_err, av_abs_err, min_dist_map


def check_frequencies(freq_truth, ab_est, truth2aln, contigs, ground_truth,
        min_HC, outfile):
    # read true frequencies
    true_frequencies = {}
    with open(freq_truth, 'r') as f:
        for line in f:
            [seq_id, freq] = line.strip('\n').split()
            true_frequencies[seq_id] = float(freq)

    total_abundance = sum(ab_est.values())

    # compute relative true frequencies, leaving out any missing haps
    total_true_freqs = 0
    for truth_id, freq in true_frequencies.items():
        aln_list = truth2aln[truth_id]
        truth_len = len(ground_truth[truth_id])
        aln_lengths = [len(contigs[aln[0][0]]) for aln in aln_list]
        if len(aln_list) > 0 and max(aln_lengths) > min_HC*truth_len:
            total_true_freqs += freq
    if total_true_freqs == 0:
        print("No contigs of sufficient length, can't evaluate frequencies.\n")
        return -1, -1, -1

    total_ab = 0
    for truth_id, aln_list in truth2aln.items():
        truth_len = len(ground_truth[truth_id])
        for info in aln_list:
            contig_id = info[0][0]
            contig_ab = ab_est[contig_id]
            if len(contigs[contig_id]) > min_HC*truth_len:
                total_ab += contig_ab
    if total_ab == 0:
        print("No contigs of sufficient length, can't evaluate frequencies.\n")
        return -1, -1, -1

    # check assignments and evaluate
    err_list = []
    abs_err_list = []
    contigs_seen = []
    f = open(outfile, 'a')
    for truth_id, aln_list in truth2aln.items():
        truth_len = len(ground_truth[truth_id])
        true_freq = true_frequencies[truth_id]
        if len(aln_list) == 0: # don't evaluate missing strains as wrong estimation
            continue
        total_ab_est = 0
        for info in aln_list:
            aln = info[0]
            contig_id = aln[0]
            if contig_id in contigs_seen:
                print("WARNING: duplicate contig assignment")
                continue
            else:
                contigs_seen.append(contig_id)
            contig_len = len(contigs[contig_id])
            contig_ab = ab_est[contig_id]
            if contig_len > min_HC*truth_len:
                # full length contig -> add estimated abundances
                total_ab_est += contig_ab
            else:
                # if not full length, evaluate frequencies manually
                print("WARNING: contig not full length")
                continue
                # return -1, -1
        freq_est = total_ab_est/total_ab*100
        cor_true_freq = true_freq/total_true_freqs*100
        if total_ab_est > 0:
            print("{}\t{}".format(cor_true_freq, freq_est))
            f.write("{}\t{}\n".format(cor_true_freq, freq_est))
            abs_err = abs(freq_est - cor_true_freq)
            rel_err = abs(freq_est - cor_true_freq)/cor_true_freq
            err_list.append(rel_err)
            abs_err_list.append(abs_err)
        #print truth_id, true_freq, freq_est, len(aln_list)
    f.close()
    if len(abs_err_list) == 1:
        print("Only 1 strain reconstructed, hence perfect frequency estimation.")
        return -1, -1, -1
    average_rel_err = sum(err_list)/len(err_list)*100
    average_abs_err = sum(abs_err_list)/len(abs_err_list)
    print("average rel error: {:.1f}%".format(average_rel_err))
    print("average abs error: {:.1f}%".format(average_abs_err))
    print("min error: {:.3f}".format(min(err_list)))
    print("max error: {:.3f}".format(max(err_list)))
    median = np.median(np.array(sorted(err_list)))
    return average_rel_err, average_abs_err, median

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
                        ab = float(line.lstrip('>').rstrip('\n').split()[-1].lstrip('frequency='))
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

def get_assembly_stats(truth2aln, total_assembly_len, N50_all_contigs,
            ground_truth, ncontigs, outfile, method):
    # prints all assembly statistics of interest
    mismatch_count = 0
    ins_count = 0
    ins_len = 0
    del_count = 0
    del_len = 0
    N_count = 0
    target_cov_len = 0 # target bases covered
    total_target_len = 0
    breakpoints = 0
    true_positives = set()
    exactly_matched = 0
    aln_contig_lengths = []
    min_dist_map = {}
    f = open(outfile, 'a')
    for truth_id, aln in truth2aln.items():
        print('\n' + truth_id + ':',)
        if len(aln) == 0:
            truth_seq_len = len(ground_truth[truth_id])
            breakpoints = 0
        else:
            truth_seq_len = aln[0][2][2]
            breakpoints += len(aln)-1
        total_target_len += truth_seq_len
        target_cov = [0 for i in range(truth_seq_len)]
        local_mismatches = 0
        local_ins_len = 0
        local_del_len = 0
        min_dist = 1
        for info in aln:
            seq_id = info[0][0]
            print("{} ({});".format(seq_id, info[0][4]),)
            scores = info[1]
            mismatch_count += scores[1]
            local_mismatches += scores[1]
            ins_count += scores[2]
            ins_len += scores[3]
            local_ins_len += scores[3]
            del_count += scores[4]
            del_len += scores[5]
            local_del_len += scores[5]
            N_count += scores[6]
            aln_range = info[2]
            aln_contig_lengths.append(aln_range[1] - aln_range[0])
            for i in range(aln_range[0], aln_range[1]):
                target_cov[i] = 1
            dist = scores[1] + scores[3] + scores[5]
            if dist == 0:
                true_positives.add(seq_id)
            assignment = [
                method, seq_id, truth_id, aln_range[1] - aln_range[0], dist,
                scores[1], scores[3], scores[5]
            ]
            f.write('\t'.join([str(x) for x in assignment]) + '\n')
            min_dist = min(min_dist, dist/(aln_range[1] - aln_range[0]))
        target_cov_len += sum(target_cov)
        min_dist_map[truth_id] = min_dist*100
        if min_dist == 0:
            exactly_matched += 1
        print("\ntarget coverage: {} of {} ({:.1f}%)".format(sum(target_cov), truth_seq_len, 100*sum(target_cov)/truth_seq_len))
        print("# contigs: {}".format(len(aln)))
        print("minimal contig distance: {:.2f}%".format(min_dist*100))
        print("mismatches, ins_len, del_len: {} {} {}".format(local_mismatches, local_ins_len, local_del_len))
    f.close()
    print()
    print("-----------------------")
    print("- Assembly statistics -")
    print("-----------------------")
    print("Total assembly length:\t{} bp".format(total_assembly_len))
    if total_assembly_len == 0:
        outfile_line = "{0}\t{0}\t{0}\t{0}\t{0}\t{0}\t{0}".format(0)
        sensitivity = 0
        ppv = 0
        return outfile_line, sensitivity, ppv, min_dist_map
    total_aln_len = sum(aln_contig_lengths)
    total_unaligned_len = total_assembly_len - total_aln_len
    unaligned_perc = float(total_unaligned_len)/total_assembly_len*100.0
    print("Total unaligned length:\t{} bp ({:.1f}%)".format(total_unaligned_len, unaligned_perc))
    edit_distance = float(mismatch_count + ins_len + del_len)/total_aln_len*100 if total_aln_len > 0 else 0
    print("Overall edit distance:\t{:5.3f}%".format(edit_distance))
    mismatch_rate = float(mismatch_count)/total_aln_len*100 if total_aln_len > 0 else 0
    print("Mismatch rate:\t{:5.3f}%".format(mismatch_rate))
    N_rate = float(N_count)/total_aln_len*100 if total_aln_len > 0 else 0
    print("'N' rate:\t{:5.3f}%".format(N_rate))
    insertion_rate = float(ins_len)/total_aln_len*100 if total_aln_len > 0 else 0
    print("Total insertion count:\t{}".format(ins_count))
    print("Total insertion length:\t{} bp ({:5.3f}%)".format(ins_len, insertion_rate))
    #print "Total insertion length:\t{:5.3f}%".format(insertion_rate)
    deletion_rate = float(del_len)/total_aln_len*100 if total_aln_len > 0 else 0
    print("Total deletion count:\t{}".format(del_count))
    print("Total deletion length:\t{} bp ({:5.3f}%)".format(del_len, deletion_rate))
    #print "Total deletion length:\t{:5.3f}%".format(deletion_rate)

    total_target_cov = float(target_cov_len)/total_target_len*100
    print("Target coverage:\t{:4.1f}%".format(total_target_cov))

    N50_aln = compute_NX(aln_contig_lengths, 50)
    print("N50 aligned sequence:\t{}".format(N50_aln))
    print("N50 all contigs:\t{}".format(N50_all_contigs))
    # TODO: minimum contig size to cover at least 50% of the ground truth

    print("Breakpoint number:\t{}".format(breakpoints))
    # TODO: number of unnecessary breakpoints
    # TODO: number of conflict cliques larger than known ploidy
    print("# true positives = {}".format(len(true_positives)))
    sensitivity = exactly_matched/len(truth2aln)
    ppv = len(true_positives)/ncontigs
    print("Sensitivity = {:.2f}".format(sensitivity))
    print("PPV = {:.2f}".format(ppv))
    print()
    outfile_line = "{}\t{}\t{}\t{}\t{:.3f}\t{}\t{:.5f}".format(N50_all_contigs,
        total_aln_len, total_target_len, target_cov_len, total_target_cov/100,
        mismatch_count+ins_len+del_len, edit_distance/100)
    return outfile_line, sensitivity, ppv, min_dist_map

def compute_NX(contig_lengths, X):
    assert X > 0
    assert X <= 100
    factor = X/100.0
    contig_lengths.sort(reverse=True)
    total_len = sum(contig_lengths)
    current_len = 0
    NX = 0
    for l in contig_lengths:
        current_len += l
        if current_len >= factor*total_len:
            NX = l
            break
    assert NX >= 0
    return NX


if __name__ == '__main__':
    sys.exit(main())
