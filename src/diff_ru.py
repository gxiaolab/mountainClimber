#!/user/bin/python -tt
"""
Calculate the relative usage of each 5' and 3' segment
"""


import os
import re
import sys
import glob
import math
import argparse
import pybedtools as pb
import numpy as np 					# v1.10.4
from scipy import stats 			# v0.15.1
from collections import *
from datetime import datetime
from functions import sort_bedfile


def get_seg2cov(infile, sample, seg2cov):
    """
    Read infile.
    Get coverage for each segment.
    """
    with open(infile, 'r') as f:
        for line in f:
            x = line.rstrip().split('\t')
            if len(x) == 6:
                (chrom, start, end, name, cov, strand) = x
            elif len(x) == 5:
                (chrom, start, end, name, cov) = x
                strand = name.split(':')[5]  # inferred strand
            gene = ':'.join(name.split(':')[1:-3])
            if (gene, chrom, start, end, strand, sample) not in seg2cov:
                seg2cov[(gene, chrom, start, end, strand, sample)] = float(cov)
            else:
                sys.stderr.write('EXIT: seen segment!: ' + ' '.join([gene, chrom, start, end, strand, sample]) + '\n')
                sys.exit(1)
    return seg2cov


def calculate_relative_usage(seg2covlist, n):
    """Calculate relative usage for each sample"""
    s2ru = {}
    for s in range(n):
        sample_cov = [seg2covlist[i][s] for i in seg2covlist]
        prxl_cov = max(sample_cov)
        if prxl_cov != 0:  # exclude samples with proximal coverage = 0 to avoid division by zero below
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


def write_output_segments(o, gene, seg, n, ind, ru):
    """Write output: segments"""
    gs, gstart, gend, chrom, strand_inferred = gene.split(':')
    (cov_mean, cov_var, label, start, end) = seg
    if strand_inferred != 'NA':
        o.write('\t'.join(map(str, [chrom, start, end, ':'.join(map(str, [label, gene, cov_mean, cov_var, n, ind])), ru, strand_inferred])) + '\n')
    else:
        o.write('\t'.join(map(str, [chrom, start, end, ':'.join(map(str, [label, gene, cov_mean, cov_var, n, ind])), ru])) + '\n')


def write_output_cp_left(o, gene, seg, n, ind, ru):
    """Write output: left end (3' for - strand; 5' for + strand)"""
    gs, gstart, gend, chrom, strand_inferred = gene.split(':')
    (cov_mean, cov_var, label, start, end) = seg
    if strand_inferred == 'NA':
        o.write('\t'.join(map(str, [chrom, int(start) - 1, start, ':'.join(map(str, [label.split('|')[0], gene, cov_mean, cov_var, n, ind])), ru])) + '\n')
    else:
        o.write('\t'.join(map(str, [chrom, int(start) - 1, start, ':'.join(map(str, [label.split('|')[0], gene, cov_mean, cov_var, n, ind])), ru, strand_inferred])) + '\n')


def write_output_cp_right(o, gene, seg, n, ind, ru):
    """Write output: right end (3' for + strand; 5' for - strand)"""
    gs, gstart, gend, chrom, strand_inferred = gene.split(':')
    (cov_mean, cov_var, label, start, end) = seg
    if strand_inferred == 'NA':
        o.write('\t'.join(map(str, [chrom, int(end) - 1, end, ':'.join(map(str, [label.split('|')[1], gene, cov_mean, cov_var, n, ind])), ru])) + '\n')
    else:
        o.write('\t'.join(map(str, [chrom, int(end) - 1, end, ':'.join(map(str, [label.split('|')[1], gene, cov_mean, cov_var, n, ind])), ru, strand_inferred])) + '\n')


def relabel_seg(seg, strand, side):
    """Relabel Junction, Intron, and Exon change points with APA or ATSS"""
    label1, label2 = seg[2].split('|')
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

    seg[2] = '|'.join(new_label_list)
    return(seg)


def main(argv):
    # --------------------------------------------------
    # get args
    # --------------------------------------------------
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Calculate the relative usage of each 5\' and 3\' segment')
    group = parser.add_argument_group('Input')
    group.add_argument('-i', '--input', dest='input', type=str, nargs='*', metavar='',
                       help='List of space-delimited output files from diff_segmentReadCounts for a single condition.')
    group.add_argument('-s', '--segments', dest='segments', type=str, metavar='',
                       help='Condition-specific _segments.bed output file from diff_cluster.')
    group.add_argument('-c', '--condition', dest='condition', type=str, metavar='',
                       help='Condition label')
    group.add_argument('-l', '--input_cp', dest='input_cp', type=str, metavar='',
                       help='Condition-specific _cluster.bed output file from diff_cluster.')

    group = parser.add_argument_group('Parameters')
    group.add_argument('-n', '--min_segments', dest='min_segments', type=int, default=3, metavar='',
                       help='Minimum number of segments required in the TU to calculate relative end usage')

    group = parser.add_argument_group('Output')
    group.add_argument('-o', '--output', dest='output', type=str, metavar='',
                       help='Output prefix.')
    group.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                       help='Print progress.')
    args = parser.parse_args()
    print args

    # --------------------------------------------------
    # main routine
    # --------------------------------------------------
    print '\njob starting:', str(datetime.now().time())

    if not args.input:
        sys.stderr.write('EXIT: Please provide --input')
        sys.exit(1)

    if not args.segments:
        sys.stderr.write('EXIT: Please provide --segments')
        sys.exit(1)

    if not args.condition:
        sys.stderr.write('EXIT: Please provide --conditions')
        sys.exit(1)

    if not args.input_cp:
        sys.stderr.write('EXIT: Please provide --input_cp')
        sys.exit(1)

    if not args.output:
        sys.stderr.write('EXIT: Please provide --output')
        sys.exit(1)

    # === I/O ===
    stdout = args.output + '_' + args.condition + '_ru_totals.txt'

    # === get genes for which we called change points ===
    genes_kept = {}
    with open(args.input_cp, 'r') as f:
        for line in f:
            x = line.rstrip().split('\t')
            gene = ':'.join(x[3].split(':')[1:6])
            genes_kept[gene] = 1

    # === get segment coverage for each sample ===
    print '- reading segment coverage', str(datetime.now().time())
    seg2cov = {}
    samples = []
    for infile in args.input:
        sample = '_'.join(os.path.basename(infile).replace(args.output + '_', '').replace('_readCounts.bed', '').split('_'))
        samples.append(sample)
        seg2cov = get_seg2cov(infile, sample, seg2cov)

    # === read segments ===
    seg_all = []
    gene2dbmax = {}
    seen = {}
    with open(args.segments, 'r') as f:
        for line in f:
            if not line.startswith('track'):
                x = line.rstrip().split('\t')
                if len(x) == 6:
                    chrom, start, end, name, score, strand = x
                elif len(x) == 5:
                    chrom, start, end, name, score = x
                    strand = name.split(':')[5]  # inferred strand
                else:
                    sys.stderr.write('EXIT: do not recognize bed file: ' + args.segments + '\n')
                    sys.exit(1)

                cplabel, gs, gstart, gend, gchrom, gstrand_inferred, winsize = name.split(':')
                cplabel_prioritized = '|'.join([x.split(';')[0] for x in cplabel.split('|')]) if ';' in cplabel else cplabel
                condLabel = '|'.join([x.split(';')[1] for x in cplabel.split('|')])
                gene = ':'.join([gs, gstart, gend, gchrom, gstrand_inferred])
                if gene in genes_kept:
                    seg_all.append((gene, chrom, start, end, gstrand_inferred, cplabel_prioritized, condLabel, score))
                    if gene not in gene2dbmax:
                        gene2dbmax[gene] = score
                    elif int(score) > int(gene2dbmax[gene]):
                        gene2dbmax[gene] = score

    # === get coverage of this segment in each condition ===
    gene2cps = defaultdict(list)
    for i, (gene, chrom, start, end, gstrand_inferred, cplabel, condLabel, ind) in enumerate(seg_all):
        ind_max = gene2dbmax[gene]
        cov = [seg2cov.get((gene, chrom, start, end, gstrand_inferred, sample), 0.0) for sample in samples]
        gene2cps[gene].append((np.mean(cov), np.var(cov), start, end, cplabel, condLabel, ind, ind_max))

    # === get relative usage for each UTR ===
    print '- calculating relative usage', str(datetime.now().time())
    o4 = open(stdout, 'w')
    o4.write('\t'.join(['Condition', 'Total genes', 'Total genes with > 3 segments', 'Total genes with change points: left end', 'Total genes with change points: right end', 'Total genes with all prxlCov != 0: left end', 'Total genes with all prxlCov != 0: right end']) + '\n')

    count_cond_gene_3segments = {}
    count_cond_gene_cps = {}
    count_cond_gene_cps['L'] = {}
    count_cond_gene_cps['R'] = {}
    count_cond_gene_filtered = {}
    count_cond_gene = {}
    count_cond_gene_min_prxlCov = {}
    count_cond_gene_min_prxlCov['L'] = {}
    count_cond_gene_min_prxlCov['R'] = {}
    outfile_seg_cond_list = []

    outfile_seg_cond = args.output + '_ru_segments_' + args.condition + '.bed'
    outfile_cp_cond = args.output + '_ru_cp_' + args.condition + '.bed'
    outfile_seg_cond_list.append(outfile_seg_cond)
    o2 = open(outfile_seg_cond, 'w')
    o3 = open(outfile_cp_cond, 'w')

    for gene in gene2cps:
        gs, gstart, gend, chrom, strand_inferred = gene.split(':')
        gene2cps_tuple = gene2cps[gene]
        cov_mean_list, cov_var_list, start_list, end_list, cplabel_list, condLabel_list, ind_list, ind_max_list = zip(*gene2cps_tuple)
        nsamples = len(samples)

        if args.verbose:
            print 'gene', gene, ind_max_list[0]
            print 'starts', start_list
            print 'ends', end_list
            print 'cplabels', cplabel_list
            print 'condLabels', condLabel_list
            print 'inds', ind_list
            print 'cov_mean', cov_mean_list

        if gene not in count_cond_gene:
            count_cond_gene[gene] = 1
        else:
            count_cond_gene[gene] += 1

        if len(cplabel_list) > args.min_segments:
            if gene not in count_cond_gene_3segments:
                count_cond_gene_3segments[gene] = 1
            else:
                count_cond_gene_3segments[gene] += 1

            # === if mountainClimber_relabel.py was used to re-annotate change point labels (e.g. to TSS_perSample or POLYA_perSample), then replace those so they don't get assigned to either end ===
            cplabel_list = list(cplabel_list)
            for c, cplabel in enumerate(cplabel_list):
                cpleft, cpright = cplabel.split('|')
                cpleft = 'INTERNAL' if '_perSample' in cpleft else cpleft
                cpright = 'INTERNAL' if '_perSample' in cpright else cpright
                cplabel_list[c] = '|'.join([cpleft, cpright])

            # === get indices for left & right ends: 1st & last exon ===
            cplabel_list_left = [x.split('|')[0] for x in cplabel_list]
            cplabel_list_right = [x.split('|')[1] for x in cplabel_list]
            if any('Junction' in x for x in cplabel_list_left) and any('Junction' in x for x in cplabel_list_right):  # multi-exon gene -> include all change points in 1st & last exon
                this_list_left = [i for i, x in enumerate(cplabel_list_left) if 'Junction' in x]
                this_list_right = [i for i, x in enumerate(cplabel_list_right) if 'Junction' in x]
                first_jxn_ind = min(this_list_left) if len(this_list_left) > 0 else 1
                last_jxn_ind = max(this_list_right) if len(this_list_right) > 0 else len(cov_mean_list) - 2
            elif any('Exon' in x and 'Tandem' not in x for x in cplabel_list):  # non-strand-specific/single exon ->  get all TANDEM points
                this_list_left = [i for i, x in enumerate(cplabel_list_left) if 'Exon' in x]
                this_list_right = [i for i, x in enumerate(cplabel_list_right) if 'Exon' in x]
                first_jxn_ind = min(this_list_left) if len(this_list_left) > 0 else 1
                last_jxn_ind = max(this_list_right) if len(this_list_right) > 0 else len(cov_mean_list) - 2
            else:  # introns present but no junctions detected --> introns are likely either short or retained. just keep first & last change point
                first_jxn_ind = 1
                last_jxn_ind = len(cov_mean_list) - 2

            indices_left = [i for i, x in enumerate(cplabel_list) if i < first_jxn_ind or ind_list[i] == '1' or 'Left' in x or ('TSS' in x and strand_inferred == '+') or ('APA' in x and strand_inferred == '-') or ('PolyA' in x and strand_inferred == '-')]
            indices_right = [i for i, x in enumerate(cplabel_list) if i > last_jxn_ind or ind_list[i] == ind_max_list[i] or 'Right' in x or ('TSS' in x and strand_inferred == '-') or ('APA' in x and strand_inferred == '+') or ('PolyA' in x and strand_inferred == '+')]

            # === check if there are any internal condition-specific change points ===
            if any('Junction' in x for x in cplabel_list):  # don't check this for single exon genes
                condLabel_list_left = [x.split('|')[0] for x in condLabel_list]
                condLabel_list_right = [x.split('|')[1] for x in condLabel_list]
                cond_indices_left = [i for i, x in enumerate(condLabel_list_left) if args.condition in x and first_jxn_ind < i < last_jxn_ind and i not in indices_left]
                cond_indices_right = [i for i, x in enumerate(condLabel_list_right) if args.condition in x and first_jxn_ind < i < last_jxn_ind and i not in indices_right]
                if args.verbose:
                    print 'cond_indices_left', cond_indices_left, [(cplabel_list_left[x], condLabel_list_left[x]) for x in cond_indices_left]
                    print 'cond_indices_right', cond_indices_right, [(cplabel_list_right[x], condLabel_list_right[x]) for x in cond_indices_right]

                if len(cond_indices_left) > 0 and cplabel_list_left[cond_indices_left[0]] != 'Junction':
                    if args.verbose:
                        print 'add change point to left!', gene, cond_indices_left[0], cplabel_list_left[cond_indices_left[0]], condLabel_list_left[cond_indices_left[0]]
                    indices_left.append(cond_indices_left[0])
                if len(cond_indices_right) > 0 and cplabel_list_right[cond_indices_right[-1]] != 'Junction':
                    if args.verbose:
                        print 'add change point to right!', gene, cond_indices_right[-1], cplabel_list_right[cond_indices_right[-1]], condLabel_list_right[cond_indices_right[-1]]
                    indices_right = [cond_indices_right[-1]] + indices_right

            # assign segments common to both ends to either end by checking whether it consecutively follows a segment at either end
            indices_common = [x for x in indices_left if x in indices_right]
            if len(indices_common) != 0:
                for x in indices_common:
                    if args.verbose:
                        print gene
                        print start_list
                        print end_list
                        print cplabel_list
                        print 'l:', indices_left
                        print 'r:', indices_right
                        print 'c:', indices_common

                    if x == 0:
                        del indices_right[indices_right.index(x)]
                        if args.verbose:
                            print 'del from right end:', x
                    elif x == len(start_list) - 1:
                        del indices_left[indices_left.index(x)]
                        if args.verbose:
                            print 'del from left end:', x
                    else:
                        if indices_left.index(x) != 0 and indices_right.index(x) != len(indices_right) - 1:
                            if indices_left[indices_left.index(x) - 1] != x - 1 and indices_right[indices_right.index(x) + 1] == x + 1:
                                if args.verbose:
                                    print 'del from left:', x
                                del indices_left[indices_left.index(x)]
                            elif indices_left[indices_left.index(x) - 1] == x - 1 and indices_right[indices_right.index(x) + 1] != x + 1:
                                if args.verbose:
                                    print 'del from right:', x
                                del indices_right[indices_right.index(x)]
                            else:
                                if indices_common.index(x) != 0 and indices_common[indices_common.index(x) - 1] != x - 1:
                                    del indices_left[indices_left.index(x)]
                                    if args.verbose:
                                        print 'del from left common:', x, indices_common[indices_common.index(x) - 1], x - 1
                                elif indices_common.index(x) != len(indices_common) - 1 and indices_common[indices_common.index(x) + 1] != x + 1:
                                    del indices_right[indices_right.index(x)]
                                    if args.verbose:
                                        print 'del from right common:', x, indices_common[indices_common.index(x) + 1], x + 1
                                else:
                                    if args.verbose:
                                        print 'del from both:', x
                                    del indices_left[indices_left.index(x)]
                                    del indices_right[indices_right.index(x)]
                        else:
                            sys.stderr.write('EXIT: can\'t assign this change point to either end\n')
                            sys.exit(1)

                if args.verbose:
                    print 'final l:', indices_left, [cplabel_list[x] for x in indices_left]
                    print 'final r:', indices_right, [cplabel_list[x] for x in indices_right]

            if args.verbose:
                print first_jxn_ind, last_jxn_ind
                print 'indices_left', indices_left, [start_list[i] for i in indices_left], [end_list[i] for i in indices_left]
                print 'indices_right', indices_right, [start_list[i] for i in indices_right], [end_list[i] for i in indices_right]
                print [condLabel_list[i] for i in indices_left]
                print [condLabel_list[i] for i in indices_right]

        else:  # did not meet minimum # segments to check for all change points -> just keep distal ends
            indices_left = [0]
            indices_right = [len(cplabel_list) - 1]

        # === distal -> proximal order ===
        segs_left_temp = []
        if len(indices_left) > 0:
            # identify alternative ends
            for i in indices_left:
                if i == 0:
                    segs_left_temp.append([cov_mean_list[i], cov_var_list[i], cplabel_list[i], start_list[i], end_list[i]])
                else:
                    if len(segs_left_temp) == 0:
                        segs_left_temp.append([cov_mean_list[i - 1], cov_var_list[i], cplabel_list[i - 1], start_list[i - 1], end_list[i - 1]])
                    segs_left_temp.append([cov_mean_list[i], cov_var_list[i], cplabel_list[i], start_list[i], end_list[i]])

        # proximal -> distal order
        segs_right_temp = []
        if len(indices_right) > 0:
            for i in indices_right:
                segs_right_temp.append([cov_mean_list[i], cov_var_list[i], cplabel_list[i], start_list[i], end_list[i]])

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

        # === remove ambiguous segments (right before proximal left & vice versa) ===
        segs_right_starts = [int(x[3]) for x in segs_right]
        segs_left_starts = [int(x[3]) for x in segs_left]

        to_del_right = []
        for i, this_start in enumerate(segs_right_starts):
            if any(x > this_start for x in segs_left_starts):
                to_del_right.append(i)
        segs_right = [x for i, x in enumerate(segs_right) if i not in to_del_right]
        # del indices_right[indices_right.index(x)]

        to_del_left = []
        for i, this_end in enumerate(segs_left_starts):
            if any(x < this_end for x in segs_right_starts):
                to_del_left.append(i)
        segs_left = [x for i, x in enumerate(segs_left) if i not in to_del_left]

        if args.verbose:
            print '-> removed ambiguous segments:'
            print 'segs_left_temp', segs_left_temp
            print 'segs_right_temp', segs_right_temp
            print segs_left
            print segs_right

        # === calculate pi and write output ===
        if len(segs_left) == 0:
            print 'NO LEFT SEG?'
            sys.exit(1)
        elif len(segs_left) == 1:
            ind = 'L0'
            this_seg = relabel_seg(segs_left[0], strand_inferred, ind)
            write_output_segments(o2, gene, this_seg, n=nsamples, ind=ind, ru=1)
            write_output_cp_left(o3, gene, this_seg, n=nsamples, ind=ind, ru=1)
        else:
            if gene not in count_cond_gene_cps['L']:
                count_cond_gene_cps['L'][gene] = 1
            else:
                count_cond_gene_cps['L'][gene] += 1

            # get coverage of each segment in each sample
            seg_left2covlist = {}
            for i, seg in enumerate(segs_left):
                (cov_mean, cov_var, label, start, end) = seg
                covlist = [seg2cov.get((gene, chrom, start, end, strand_inferred, sample), 0.0) for sample in samples]
                seg_left2covlist[i] = covlist
                if i > 0:
                    if np.mean(seg_left2covlist[i]) < np.mean(seg_left2covlist[i - 1]):
                        covlist = [0] * len(samples)
                        seg_left2covlist[i] = covlist

            # no proximal expression = 0
            if gene not in count_cond_gene_min_prxlCov['L']:
                count_cond_gene_min_prxlCov['L'][gene] = 1
            else:
                count_cond_gene_min_prxlCov['L'][gene] += 1

            ru = calculate_relative_usage(seg_left2covlist, len(samples))  # , segs_left[-1][0])
            if args.verbose:
                print 'ru left', ru

            segs_left_nonzeroRU = [x for i, x in enumerate(segs_left) if ru[i] != 0]
            nonzeroRU = [x for x in ru if x != 0]
            for i, x in enumerate(segs_left_nonzeroRU):
                ind = 'L' + str(i)
                this_seg = relabel_seg(x, strand_inferred, ind)
                write_output_segments(o2, gene, this_seg, n=nsamples, ind=ind, ru=nonzeroRU[i])
                write_output_cp_left(o3, gene, this_seg, n=nsamples, ind=ind, ru=nonzeroRU[i])

        if len(segs_right) == 0:
            print 'NO RIGHT SEG?'
            sys.exit(1)
        elif len(segs_right) == 1:
            ind = 'R0'
            this_seg = relabel_seg(segs_right[-1], strand_inferred, ind)
            write_output_segments(o2, gene, this_seg, n=nsamples, ind=ind, ru=1)
            write_output_cp_right(o3, gene, this_seg, n=nsamples, ind=ind, ru=1)
        else:
            if gene not in count_cond_gene_cps['R']:
                count_cond_gene_cps['R'][gene] = 1
            else:
                count_cond_gene_cps['R'][gene] += 1

            # get coverage of each segment in each sample
            seg_right2covlist = {}
            for i, seg in enumerate(segs_right[::-1]):  # reverse: distal -> proximal order
                (cov_mean, cov_var, label, start, end) = seg
                covlist = [seg2cov.get((gene, chrom, start, end, strand_inferred, sample), 0.0) for sample in samples]
                seg_right2covlist[i] = covlist
                if i > 0:
                    # if mean(proximal) < mean(distal) across samples, report proximal as 0. do this instead of filtering out so that if we have alternative last exons in 2 conditions then both change points will be reported
                    if np.mean(seg_right2covlist[i]) < np.mean(seg_right2covlist[i - 1]):
                        covlist = [0] * len(samples)
                        seg_right2covlist[i] = covlist

            # no proximal expression = 0
            if gene not in count_cond_gene_min_prxlCov['R']:
                count_cond_gene_min_prxlCov['R'][gene] = 1
            else:
                count_cond_gene_min_prxlCov['R'][gene] += 1

            ru = calculate_relative_usage(seg_right2covlist, len(samples))
            ru = ru[::-1]  # reverse: proximal -> distal order
            if args.verbose:
                print 'ru right', ru

            segs_right_nonzeroRU = [x for i, x in enumerate(segs_right) if ru[i] != 0]
            nonzeroRU = [x for x in ru if x != 0]
            for i, x in enumerate(segs_right_nonzeroRU):
                ind = 'R' + str(len(segs_right_nonzeroRU) - 1 - i)
                this_seg = relabel_seg(x, strand_inferred, ind)
                write_output_segments(o2, gene, this_seg, n=nsamples, ind=ind, ru=nonzeroRU[i])
                write_output_cp_right(o3, gene, this_seg, n=nsamples, ind=ind, ru=nonzeroRU[i])
    o2.close()
    o3.close()

    o4.write('\t'.join(map(str, [args.condition, len(count_cond_gene.keys()), len(count_cond_gene_3segments.keys()), len(count_cond_gene_cps['L'].keys()), len(count_cond_gene_cps['R'].keys()), len(count_cond_gene_min_prxlCov['L'].keys()), len(count_cond_gene_min_prxlCov['R'].keys())])) + '\n')

    # sort bed files
    for file in [outfile_seg_cond, outfile_cp_cond]:
        sort_bedfile(infile=file, outfile=file + '.sorted')
        os.rename(file + '.sorted', file)
    o4.close()

    # remove temporary directory
    if os.path.isdir(args.output):
        os.rmdir(args.output)

    print 'finished:', str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
    main(sys.argv[1:])
