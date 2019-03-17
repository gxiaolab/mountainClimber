#!/user/bin/python -tt

import os
import re
import sys
import glob
import math
import argparse
import numpy as np  # v1.10.4
from scipy import stats  # v0.15.1
from collections import *
from datetime import datetime
from functions import run_command, sort_bedfile


def calculate_relative_usage(seg2covlist, n):
    """ calculate relative usage for each sample in distal to proximal order """
    s2ru = {}
    for s in range(n):
        sample_cov = [seg2covlist[i][s] for i in seg2covlist]
        prxl_cov = max(sample_cov)
        sample_cov_sorted = sorted(sample_cov)
        sample_cov_sorted_indices = [sample_cov.index(x) for x in sample_cov_sorted]
        ru_unsorted = [x / prxl_cov if j == 0 else (x - sample_cov_sorted[j - 1]) / prxl_cov for j, x in enumerate(sample_cov_sorted)]
        ru = [x for _, x in sorted(zip(sample_cov_sorted_indices, ru_unsorted))]  # sort ru by the original change point order
        s2ru[s] = ru

    # calculate mean relative usage across samples
    ru = []
    for i in seg2covlist:
        seg_ru_list = [s2ru[s][i] for s in s2ru]
        ru.append(np.mean(seg_ru_list))

    return ru


def write_output_cp_left(o, gene, seg, n, ind, ru):
    """write output: left side"""
    gs, gstart, gend, strand_inferred, chrom = gene.split(':')
    (cov_mean, label, start, end) = seg
    if strand_inferred == 'NA':
        o.write('\t'.join(map(str, [chrom, int(start) - 1, start, ':'.join(map(str, [label.split('|')[0], gene, cov_mean, n, ind])), ru])) + '\n')
    else:
        o.write('\t'.join(map(str, [chrom, int(start) - 1, start, ':'.join(map(str, [label.split('|')[0], gene, cov_mean, n, ind])), ru, strand_inferred])) + '\n')


def write_output_cp_right(o, gene, seg, n, ind, ru):
    """write output: right side"""
    gs, gstart, gend, strand_inferred, chrom = gene.split(':')
    (cov_mean, label, start, end) = seg
    if strand_inferred == 'NA':
        o.write('\t'.join(map(str, [chrom, int(end) - 1, end, ':'.join(map(str, [label.split('|')[1], gene, cov_mean, n, ind])), ru])) + '\n')
    else:
        o.write('\t'.join(map(str, [chrom, int(end) - 1, end, ':'.join(map(str, [label.split('|')[1], gene, cov_mean, n, ind])), ru, strand_inferred])) + '\n')


def relabel_seg(seg, strand, side):
    """Relabel Junction, Intron, and Exon change points with APA or ATSS"""
    label1, label2 = seg[1].split('|')
    new_label_list = []
    for this_label in [label1, label2]:
        if this_label == 'Junction':
            if strand == '+':
                new_label = 'JunctionAPA' if 'R' in side else 'JunctionATSS'
            elif strand == '-':
                new_label = 'JunctionAPA' if 'L' in side else 'JunctionATSS'
            else:
                new_label = this_label
        elif this_label == 'Intron':
            if strand == '+':
                new_label = 'IntronicAPA' if 'R' in side else 'IntronicATSS'
            elif strand == '-':
                new_label = 'IntronicAPA' if 'L' in side else 'IntronicATSS'
            else:
                new_label = this_label
        elif this_label == 'Exon':
            if strand == '+':
                new_label = 'ExonicAPA' if 'R' in side else 'ExonicATSS'
            elif strand == '-':
                new_label = 'ExonicAPA' if 'L' in side else 'ExonicATSS'
            else:
                new_label = this_label
        else:
            new_label = this_label
        new_label_list.append(new_label)

    seg[1] = '|'.join(new_label_list)
    return(seg)


def main(argv):
    # --------------------------------------------------
    # get args
    # --------------------------------------------------
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Calculate relative usage of each end isoform. Run this per condition. ')
    group = parser.add_argument_group('input')
    group.add_argument('-i', '--input', dest='input', type=str, metavar='', help='Bed file of change points.')
    group = parser.add_argument_group('Parameters')
    group.add_argument('-n', '--min_segments', dest='min_segments', type=int, default=3, metavar='', help='Minimum number of segments required in the TU to calculate relative end usage.')
    group = parser.add_argument_group('output')
    group.add_argument('-o', '--output', dest='output', type=str, help='Output bed filename. Bed name field = CPlabel:gene:TUstart:TUend:inferred_strand:chromosome:segmentCoverage:CPindex')
    group.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Print progress')
    args = parser.parse_args()
    print args

    # --------------------------------------------------
    # main routine
    # --------------------------------------------------
    print '\njob starting:', str(datetime.now().time())

    # === get coverage for each end segment ===
    seg2cov = {}
    gene2cps = defaultdict(list)
    gene2segs = defaultdict(list)
    sample = os.path.basename(args.input)
    samples = [sample]
    with open(args.input, 'r') as f:
        for line in f:
            if not line.startswith('track'):
                x = line.rstrip().split('\t')
                chrom = x[0]
                cpstart = x[1]
                cpend = x[2]
                name = x[3]
                gene = ':'.join(name.split(':')[1:5] + [chrom])
                gene2segs[gene].append((chrom, cpstart, cpend, name))

    for gene in gene2segs:
        gs, gstart, gend, strand_inferred, chrom = gene.split(':')
        for i, (chrom, cpstart, cpend, name) in enumerate(gene2segs[gene]):
            name_list = name.split(':')
            cplabel = name_list[0]
            ind_max = len(gene2segs[gene]) - 1
            condition = name_list[-1] if len(args.input) > 1 else 'NA'

            if i > 0:
                cov_mean = prev_cov_mean
                seg_start = prev_cpstart if i == 1 else prev_cpend
                seg_end = cpend
                gene2cps[gene].append((float(cov_mean), int(seg_start), int(seg_end), '|'.join([prev_cplabel, cplabel]), i, ind_max))
                seg2cov[(gene, chrom, int(seg_start), int(seg_end), strand_inferred, sample)] = float(cov_mean)

            # previous segment info
            prev_cov_mean = name_list[6]
            prev_cpstart = cpstart
            prev_cpend = cpend
            prev_cplabel = cplabel
            if len(args.input) > 1:
                prev_condition = condition

    o3 = open(args.output, 'w')

    for gene in gene2cps:
        gs, gstart, gend, strand_inferred, chrom = gene.split(':')
        gene2cps_tuple = gene2cps[gene]
        cov_mean_list, start_list, end_list, cplabel_list, ind_list, ind_max_list = zip(*gene2cps_tuple)
        nsamples = len(args.input)
        if args.verbose:
            print 'gene', gene, ind_max_list[0]
            print 'cov', cov_mean_list
            print 'starts', start_list
            print 'ends', end_list
            print 'cplabels', cplabel_list
            print 'inds', ind_list

        # if len(cplabel_list) >= 3:
        if len(cplabel_list) > args.min_segments:
            # === get indices for left & right ends: 1st & last exon ===
            cplabel_list_left = [x.split('|')[0] for x in cplabel_list]
            cplabel_list_right = [x.split('|')[1] for x in cplabel_list]
            if any('Junction' in x for x in cplabel_list_left) and any('Junction' in x for x in cplabel_list_right):  # multi-exon gene -> include all change points in 1st & last exon
                first_jxn_ind = min([i for i, x in enumerate(cplabel_list_left) if 'Junction' in x])
                last_jxn_ind = max([i for i, x in enumerate(cplabel_list_right) if 'Junction' in x])
            elif any('Exon' in x and 'Tandem' not in x for x in cplabel_list):  # non-strand-specific/single exon ->  get all TANDEM points
                first_jxn_ind = min([i for i, x in enumerate(cplabel_list_left) if 'Exon' in x])
                last_jxn_ind = max([i for i, x in enumerate(cplabel_list_right) if 'Exon' in x])
            else:  # introns present but no junctions detected --> introns are likely either short or retained. just keep first & last change point
                first_jxn_ind = 1
                last_jxn_ind = len(cov_mean_list) - 2
            if args.verbose:
                print 'jxn inds:', first_jxn_ind, last_jxn_ind

            indices_left = [i for i, x in enumerate(cplabel_list) if i < first_jxn_ind or ind_list[i] == '1' or 'Left' in x or ('TSS' in x and strand_inferred == '+') or ('APA' in x and strand_inferred == '-') or ('PolyA' in x and strand_inferred == '-')]
            indices_right = [i for i, x in enumerate(cplabel_list) if i > last_jxn_ind or ind_list[i] == ind_max_list[i] or 'Right' in x or ('TSS' in x and strand_inferred == '-') or ('APA' in x and strand_inferred == '+') or ('PolyA' in x and strand_inferred == '+')]
        else:
            indices_left = [0]
            indices_right = [len(cplabel_list) - 1]

        # === distal -> proximal order ===
        segs_left_temp = []
        if len(indices_left) > 0:
            # identify alternative ends
            for i in indices_left:
                if i == 0:
                    segs_left_temp.append([cov_mean_list[i], cplabel_list[i], start_list[i], end_list[i]])
                else:
                    if len(segs_left_temp) == 0:
                        segs_left_temp.append([cov_mean_list[i - 1], cplabel_list[i - 1], start_list[i - 1], end_list[i - 1]])
                    segs_left_temp.append([cov_mean_list[i], cplabel_list[i], start_list[i], end_list[i]])

        # proximal -> distal order
        segs_right_temp = []
        if len(indices_right) > 0:
            # identify alternative ends
            for i in indices_right:
                segs_right_temp.append([cov_mean_list[i], cplabel_list[i], start_list[i], end_list[i]])

        # === remove ambiguous segments (in both ends) ===
        segs_left = list(segs_left_temp)
        segs_right = list(segs_right_temp)
        segs_common = [x for x in segs_left if x in segs_right]
        if len(segs_common) > 1:
            for seg in segs_common:
                if seg != segs_left[0] and seg != segs_right[-1]:
                    segs_left.remove(seg)
                    segs_right.remove(seg)

        if len(segs_left) != len(segs_left_temp) or len(segs_right) != len(segs_right_temp):
            sys.stderr.write('WARNING: removed ambiguous segments that could have been assigned to either end\n')

        if args.verbose:
            print 'segs_left_temp', segs_left_temp
            print 'segs_right_temp', segs_right_temp
            print 'segs_left', segs_left
            print 'segs_right', segs_right

        # === calculate pi and write output ===
        if len(segs_left) == 0:
            print 'NO LEFT SEG?'
            sys.exit(1)
        elif len(segs_left) == 1:
            ind = 'L0'
            this_seg = relabel_seg(segs_left[0], strand_inferred, ind)
            write_output_cp_left(o3, gene, this_seg, n=nsamples, ind=ind, ru=1)
        else:
            # get coverage of each segment in each sample
            seg_left2covlist = {}
            for i, seg in enumerate(segs_left):
                (cov_mean, label, start, end) = seg
                covlist = [seg2cov.get((gene, chrom, start, end, strand_inferred, sample), 0.0) for sample in samples]
                seg_left2covlist[i] = covlist
                if i > 0:
                    if np.mean(seg_left2covlist[i - 1]) == 0:
                        covlist = [0] * len(samples)
                        seg_left2covlist[i] = covlist
                    elif np.mean(seg_left2covlist[i]) < np.mean(seg_left2covlist[i - 1]):
                        covlist = [0] * len(samples)
                        seg_left2covlist[i] = covlist

            ru = calculate_relative_usage(seg_left2covlist, len(samples))
            if args.verbose:
                print 'ru left', ru

            segs_left_nonzeroRU = [x for i, x in enumerate(segs_left) if ru[i] != 0]
            nonzeroRU = [x for x in ru if x != 0]
            for i, x in enumerate(segs_left_nonzeroRU):
                ind = 'L' + str(i)
                this_seg = relabel_seg(x, strand_inferred, ind)
                write_output_cp_left(o3, gene, this_seg, n=nsamples, ind=ind, ru=nonzeroRU[i])

        if len(segs_right) == 0:
            print 'NO RIGHT SEG?'
            sys.exit(1)
        elif len(segs_right) == 1:
            ind = 'R0'
            this_seg = relabel_seg(segs_right[-1], strand_inferred, ind)
            write_output_cp_right(o3, gene, this_seg, n=nsamples, ind=ind, ru=1)
        else:
            # get coverage of each segment in each sample
            seg_right2covlist = {}
            for i, seg in enumerate(segs_right[::-1]):  # reverse: distal -> proximal order
                (cov_mean, label, start, end) = seg
                covlist = [seg2cov.get((gene, chrom, start, end, strand_inferred, sample), 0.0) for sample in samples]
                seg_right2covlist[i] = covlist
                if i > 0:
                    if np.mean(seg_right2covlist[i - 1]) == 0:
                        covlist = [0] * len(samples)
                        seg_right2covlist[i] = covlist
                    elif np.mean(seg_right2covlist[i]) < np.mean(seg_right2covlist[i - 1]):
                        covlist = [0] * len(samples)
                        seg_right2covlist[i] = covlist

            ru = calculate_relative_usage(seg_right2covlist, len(samples))
            ru = ru[::-1]  # reverse: proximal -> distal order
            if args.verbose:
                print 'ru right', ru

            segs_right_nonzeroRU = [x for i, x in enumerate(segs_right) if ru[i] != 0]
            nonzeroRU = [x for x in ru if x != 0]
            for i, x in enumerate(segs_right_nonzeroRU):
                ind = 'R' + str(len(segs_right_nonzeroRU) - 1 - i)
                this_seg = relabel_seg(x, strand_inferred, ind)
                write_output_cp_right(o3, gene, this_seg, n=nsamples, ind=ind, ru=nonzeroRU[i])
    o3.close()

    sort_bedfile(infile=args.output, outfile=args.output + '.sorted')
    os.rename(args.output + '.sorted', args.output)

    print 'finished:', str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
    main(sys.argv[1:])
