#!/user/bin/python -tt
"""
Test for differential ATSS and APA sites between two conditions
"""


import os
import sys
import argparse
import numpy as np 					# v1.10.4
from collections import defaultdict
from datetime import datetime
import functions
from functions import run_command
import pybedtools as pb


def get_seg2cov(infile, sample, seg2cov):
	"""Get segment coverage"""
	with open(infile, 'r') as f:
		for line in f:
			x = line.rstrip().split('\t')
			if len(x) == 6:
				(chrom, start, end, name, cov, strand) = x
			elif len(x) == 5:
				(chrom, start, end, name, cov) = x
				strand = name.split(':')[5] # inferred strand
			gene = ':'.join(name.split(':')[1:-3])
			if (gene, chrom, start, end, strand, sample) not in seg2cov:
				seg2cov[(gene, chrom, start, end, strand, sample)] = float(cov)
			else:
				sys.stderr.write('EXIT: seen segment!: ' + ' '.join([gene, chrom, start, end, strand, sample]) + '\n')
				sys.exit(1)
	return seg2cov


def read_input(infile):
	"""get all change points per gene"""
	gene2seg_left = defaultdict(list)
	gene2seg_right = defaultdict(list)
	with open(infile, 'r') as f:
		for line in f:
			if not line.startswith('track'):
				x = line.rstrip().split('\t')
				chrom, seg_start, seg_end, name, ru = x[:5]
				label, gs, gstart, gend, gchrom, gstrand_inferred, cov_mean, cov_var, nsample, ind = name.split(':')
				ind_side = ind[0]
				ind_num = int(ind[1:])
				seg_start, seg_end = map(int, [seg_start, seg_end])
				cov_mean, ru = map(float, [cov_mean, ru])
				gene = ':'.join([gs, gstart, gend, chrom, gstrand_inferred])
				if ind_side == 'L':
					gene2seg_left[gene].append((seg_start, seg_end, label, cov_mean, ru, ind))
				elif ind_side == 'R':
					gene2seg_right[gene].append((seg_start, seg_end, label, cov_mean, ru, ind))
				else:
					sys.stderr.write('EXIT: change point must be labeled L or R!\n')
					sys.exit(1)

	cp2prxl = {}
	cp2dstl = {}
	gene2prxl = {}
	gene2dstl = {}

	# get proximal & distal: left = proximal before distal
	for gene in gene2seg_left:
		cstart_list, cend_list, label_list, cov_mean_list, ru_list, ind_list = zip(*sorted(gene2seg_left[gene]))
		if len(label_list) == 1:
			label = label_list[0].split('|')[0]
			cp_start = cstart_list[0] - 1
			cp_end = cstart_list[0]
			cp2dstl[(cp_start, cp_end, label, gene)] = (ind_list[0], cstart_list[0], cend_list[0],
				cov_mean_list[0], ru_list[0], label)
			gene2dstl[gene] = 1
		else:
			for i, ind in enumerate(label_list):
				if i > 0:
					cp_start = cstart_list[i] - 1
					cp_end = cstart_list[i]
					label = label_list[i].split('|')[0]
					label_dstl = label_list[i-1].split('|')[0]
					cp2dstl[(cp_start, cp_end, label, gene)] = (ind_list[i-1], cstart_list[i-1],
						cend_list[i-1], cov_mean_list[i-1], ru_list[i-1], label_dstl)
					gene2dstl[gene] = 1
					cp2prxl[(cp_start, cp_end, label, gene)] = (ind_list[i], cstart_list[i],
						cend_list[i], cov_mean_list[i], ru_list[i], label)
					gene2prxl[gene] = 1

	# get proximal & distal: left = proximal before distal
	for gene in gene2seg_right:
		cstart_list, cend_list, label_list, cov_mean_list, ru_list, ind_list = zip(*gene2seg_right[gene])
		if len(label_list) == 1:
			cp_start = cend_list[0] - 1
			cp_end = cend_list[0]
			label = label_list[0].split('|')[1]
			cp2dstl[(cp_start, cp_end, label, gene)] = (ind_list[0], cstart_list[0], cend_list[0],
				cov_mean_list[0], ru_list[0], label)
			gene2dstl[gene] = 1
		else:
			for i, ind in enumerate(label_list):
				if i < len(label_list) - 1:
					cp_start = cend_list[i] - 1
					cp_end = cend_list[i]
					label = label_list[i].split('|')[1]
					label_dstl = label_list[i+1].split('|')[1]
					cp2dstl[(cp_start, cp_end, label, gene)] = (ind_list[i+1], cstart_list[i+1],
						cend_list[i+1], cov_mean_list[i+1], ru_list[i+1], label_dstl)
					gene2dstl[gene] = 1
					cp2prxl[(cp_start, cp_end, label, gene)] = (ind_list[i], cstart_list[i],
						cend_list[i], cov_mean_list[i], ru_list[i], label)
					gene2prxl[gene] = 1

	return cp2prxl, cp2dstl, gene2dstl, gene2prxl


def test_diffl_cp(cov_prxla_list, cov_dstla_list, cov_prxlb_list, cov_dstlb_list, total_cp_min_dstlCov,
	total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing,
	total_cp_strictly_decreasing_filtered):
	"""get values for differential usage test"""
	cov_mean_dstla_scaled = 'NA'
	cov_mean_dstlb_scaled = 'NA'
	dstl_prxl_a = 'NA'
	dstl_prxl_b = 'NA'
	cv2a = 'NA'
	cv2b = 'NA'

	total_cp_min_dstlCov += 1
	max_prxl = max(cov_prxla_list + cov_prxlb_list)
	dstl_prxl_a = np.mean([(x + 1) / (cov_prxla_list[i] + 1) for i,x in enumerate(cov_dstla_list)])
	if dstl_prxl_a > 1:
		dstl_prxl_a = 'NA'
	dstl_prxl_b = np.mean([(x + 1) / (cov_prxlb_list[i] + 1) for i,x in enumerate(cov_dstlb_list)])
	if dstl_prxl_b > 1:
		dstl_prxl_b = 'NA'
	total_cp_nonzero_prxl += 1

	# require strictly decreasing from proximal to distal
	if dstl_prxl_a == 'NA' or dstl_prxl_b == 'NA':
		total_cp_strictly_decreasing_filtered += 1
	else:
		total_cp_strictly_decreasing += 1

		# scale proximal segment cov in all samples to the sample with max proximal coverage
		cov_dstla_list_scaled = [x * (max_prxl + 1) / (cov_prxla_list[j] + 1) for j,x in enumerate(cov_dstla_list)]
		cov_dstlb_list_scaled = [x * (max_prxl + 1) / (cov_prxlb_list[j] + 1) for j,x in enumerate(cov_dstlb_list)]
		cov_mean_dstla_scaled = np.mean(cov_dstla_list_scaled)
		cov_mean_dstlb_scaled = np.mean(cov_dstlb_list_scaled)
		if cov_mean_dstla_scaled != 0:
			cv2a = np.var(cov_dstla_list_scaled) / cov_mean_dstla_scaled ** 2
		if cov_mean_dstlb_scaled != 0:
			cv2b = np.var(cov_dstlb_list_scaled) / cov_mean_dstlb_scaled ** 2

	return cov_mean_dstla_scaled, cov_mean_dstlb_scaled, dstl_prxl_a, dstl_prxl_b, cv2a, cv2b, total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered


def get_test_type(strand, side, this_type):
	"""Decide whether the change point is: AFE, ALE, AE, TandemATSS, TandemAPA"""
	if this_type == 'Tandem':
		if (strand == '+' and side == 'L') or (strand == '-' and side == 'R'):
			test_type = 'TandemATSS'
		elif (strand == '+' and side == 'R') or (strand == '-' and side == 'L'):
			test_type = 'TandemAPA'
		else:
			test_type = 'Tandem'
	elif this_type == 'AE':
		if (strand == '+' and side == 'L') or (strand == '-' and side == 'R'):
			test_type = 'AFE'
		elif (strand == '+' and side == 'R') or (strand == '-' and side == 'L'):
			test_type = 'ALE'
		else:
			test_type = 'AE'
	else:
		sys.stderr.write('EXIT: do not recognize this type: ' + this_type + '\n')
		sys.exit(1)
	return(test_type)


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Test for differential ATSS and APA sites between two conditions.')
	group = parser.add_argument_group('Input')
	group.add_argument('-i', '--input', dest='input', type=str, nargs='*', metavar='',
		help='List of space-delimited output files from diff_segmentReadCounts for two conditions.')
	group.add_argument('-r', '--ru_segments', dest='ru_segments', type=str, nargs='*', metavar='',
		help='List of space-delimited _ru_segments output files from diff_ru, \
		one for each condition')
	group.add_argument('-ci', '--conditions_input', dest='conditions_input', type=str, nargs='*', metavar='',
		help='List of space-delimited condition labels for each --input file.')
	group.add_argument('-cr', '--conditions_ru_segments', dest='conditions_ru_segments', type=str, nargs='*', metavar='',
		help='List of space-delimited condition labels for each --ru_segments file.')

	group = parser.add_argument_group('Parameters')
	group.add_argument('-d', '--min_dstlCov', dest='min_dstlCov', type=float, default=5.0, metavar='',
		help='Minimum average reads per bp in distal segment across samples in at least 1 condition')
	group.add_argument('-p', '--min_prxlCov', dest='min_prxlCov', type=float, default=0.0, metavar='',
		help='Minimum average reads per bp in proximal segment across samples in all conditions')
	group.add_argument('-t', '--pmax', dest='pmax', type=float, default=0.05, metavar='',
		help='Maximum p-value.')
	group.add_argument('-m', '--dtop_abs_dif_min', dest='dtop_abs_dif_min', type=float, default=0.05, metavar='',
		help='Minimum relative usage (RU) difference.')

	group = parser.add_argument_group('Output')
	group.add_argument('-o', '--output', dest='output', type=str, metavar='',
		help='Output prefix. 3 output files: (1) bed file of all tested change points, (2) bed file of differential change points, (3) summary of total differential change points. Bed name field = CPlabel;testType;CPindex;gene:TUstart:TUend:proximalSegment:distalSegment:RUdifference:RUconditionA:RUconditionB:pBH')
	group.add_argument('-k', '--keep', dest='keep', action='store_true',
		help='Keep intermediate output files.')
	group.add_argument('-v', '--verbose', dest='verbose', action='store_true',
		help='Print progress.')
	args = parser.parse_args()
	print args

	# --------------------------------------------------
	# main routine
	# --------------------------------------------------
	print '\njob starting:', str(datetime.now().time())

	if not args.output:
		sys.stderr.write('EXIT: Please provide --output')
		sys.exit(1)

	if args.input:
		if len(args.conditions_input) != 1 and len(args.input) != len(args.conditions_input):
			sys.stderr.write(' '.join(map(str, ['WARNING: number of samples don\'t match!:', len(args.input), len(args.conditions_input)])) + '\n')
		if len(list(set(args.conditions_input))) != len(list(set(args.conditions_ru_segments))):
			sys.stderr.write('EXIT: the number of unique --conditions_input should equal the total unique --conditions_ru_segments input\n')
			sys.exit(1)

	if not args.ru_segments:
		sys.stderr.write('EXIT: Please provide --ru_segments')
		sys.exit(1)
	elif len(args.ru_segments) != 1:
		if len(args.ru_segments) != len(list(set(args.conditions_ru_segments))):
			sys.stderr.write('EXIT: the number of --ru_segments files should equal the total unique --conditions input\n')
			sys.exit(1)

	# === get segment coverage for each sample ===
	if len(args.ru_segments) >= 2: # if only one RU segment input, don't need to re-calculate read counts for the segments. in this case, we just label tandem vs. alternative first/last exon
		print '- reading segment coverage', str(datetime.now().time())
		seg2cov = {}
		samples = []
		cond2samples = defaultdict(list)
		for c, cond in enumerate(args.conditions_input):
			sample = os.path.basename(args.input[c]).replace('_readCounts.bed', '')
			samples.append(sample)
			cond2samples[cond].append(sample)
			print args.input[c], sample
			seg2cov = get_seg2cov(args.input[c], sample, seg2cov)

	# --------------------------------------------------
	# two conditions
	# --------------------------------------------------
	if len(args.ru_segments) == 2:
		# === I/O ===
		outfile_cp_allTested = args.output + '_cp_allTested.bed'
		outfile_cp_diff = args.output + '_cp_diff.bed'
		file_test = args.output + '_test.txt'
		file_totals = args.output + '_test_totals.txt'

		# === read input ===
		cp2prxla, cp2dstla, gene2dstla, gene2prxla = read_input(args.ru_segments[0])
		cp2prxlb, cp2dstlb, gene2dstlb, gene2prxlb = read_input(args.ru_segments[1])

		# === test before vs. after each change point ===
		o5 = open(file_test, 'w')
		o5.write('cp\tgene\ttest_type\tside\tseg_prxl\tseg_dstl\tru_dif\tru_prxla\tru_prxlb\tru_dstla\tru_dstlb\tmean_a\tmean_b\tcv2_a\tcv2_b\tnsamples_a\tnsamples_b\n')

		samplesa = cond2samples[cond2samples.keys()[0]]
		samplesb = cond2samples[cond2samples.keys()[1]]

		tested_count = 0
		total_cp = 0
		total_cp_min_dstlCov = 0
		total_cp_min_dstlCov_filtered = 0
		total_cp_nonzero_prxl_filtered = 0
		total_cp_nonzero_prxl = 0
		total_cp_strictly_decreasing = 0
		total_cp_strictly_decreasing_filtered = 0

		gene2diffcp = defaultdict(list)
		for (cp_start, cp_end, label, gene) in sorted(list(set(cp2dstla.keys() + cp2dstlb.keys()))):
			total_cp += 1
			chrom = gene.split(':')[3]
			strand = gene.split(':')[-1]

			if args.verbose:
				print cp_start, cp_end, label, gene

			# === get prxl/dstl info and coverage in each condition ===
			if (cp_start, cp_end, label, gene) in cp2prxla:
				side_prxla, cstart_prxla, cend_prxla, cov_mean_prxla, ru_prxla, label_prxla = cp2prxla[(cp_start, cp_end, label, gene)]
				if args.verbose:
					print 'prxla', side_prxla, cstart_prxla, cend_prxla, cov_mean_prxla, ru_prxla
			if (cp_start, cp_end, label, gene) in cp2prxlb:
				side_prxlb, cstart_prxlb, cend_prxlb, cov_mean_prxlb, ru_prxlb, label_prxlb = cp2prxlb[(cp_start, cp_end, label, gene)]
				if args.verbose:
					print 'prxlb', side_prxlb, cstart_prxlb, cend_prxlb, cov_mean_prxlb, ru_prxlb
			if (cp_start, cp_end, label, gene) in cp2dstla:
				side_dstla, cstart_dstla, cend_dstla, cov_mean_dstla, ru_dstla, label_dstla = cp2dstla[(cp_start, cp_end, label, gene)]
				if args.verbose:
					print 'dstla', side_dstla, cstart_dstla, cend_dstla, cov_mean_dstla, ru_dstla
			if (cp_start, cp_end, label, gene) in cp2dstlb:
				side_dstlb, cstart_dstlb, cend_dstlb, cov_mean_dstlb, ru_dstlb, label_dstlb = cp2dstlb[(cp_start, cp_end, label, gene)]
				if args.verbose:
					print 'dstlb', side_dstlb, cstart_dstlb, cend_dstlb, cov_mean_dstlb, ru_dstlb

			# === get differential change points ===
			if (cp_start, cp_end, label, gene) in cp2dstla:
				if (cp_start, cp_end, label, gene) in cp2prxla:
					if (cp_start, cp_end, label, gene) in cp2dstlb:
						if (cp_start, cp_end, label, gene) in cp2prxlb:

							# === quality check ===
							if ((cp_start, cp_end, label, gene) in cp2prxla and (cp_start, cp_end, label, gene) in cp2prxlb):
								if cstart_prxla != cstart_prxlb or cend_prxla != cend_prxlb:
									sys.stderr.write(' '.join(map(str, ['WARNING: prxl segments do not match across conditions! skipping this change point:', cp_start, cp_end, label, gene, cstart_prxla, cstart_prxlb, cend_prxla, cend_prxlb])) + '\n')
									continue
							elif ((cp_start, cp_end, label, gene) in cp2dstla and (cp_start, cp_end, label, gene) in cp2dstlb):
								if cstart_dstla != cstart_dstlb or cend_dstla != cend_dstlb:
									sys.stderr.write(' '.join(map(str, ['WARNING: dstl segments do not match across conditions! skipping this change point:', cp_start, cp_end, label, gene, cstart_dstla, cstart_dstlb, cend_dstla, cend_dstlb])) + '\n')
									continue

							if (side_prxla[0] == 'L' and cstart_prxla != cend_dstla) or (side_prxla[0] == 'R' and cend_prxla != cstart_dstla) or (side_prxlb[0] == 'L' and cstart_prxlb != cend_dstlb) or (side_prxlb[0] == 'R' and cend_prxlb != cstart_dstlb):
								if args.verbose:
									print '-> change point in both A & B: not consecutive'

								test_type = get_test_type(strand, side_prxla[0], 'AE')
								gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
									test_type, side_prxla,
									str(cstart_prxla) + '-' + str(cend_prxla),
									str(cstart_dstla) + '-' + str(cend_dstla),
									float(ru_prxla) - float(ru_prxlb),
									ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))

							else:
								# === test ===
								if args.verbose:
									print '-> change point in both A & B: test', cov_mean_dstla, cov_mean_dstlb
								cov_prxla_list = [seg2cov.get((gene, chrom, str(cstart_prxla), str(cend_prxla), strand, sample), 0.0) for sample in samplesa]
								cov_prxlb_list = [seg2cov.get((gene, chrom, str(cstart_prxlb), str(cend_prxlb), strand, sample), 0.0) for sample in samplesb]
								cov_dstla_list = [seg2cov.get((gene, chrom, str(cstart_dstla), str(cend_dstla), strand, sample), 0.0) for sample in samplesa]
								cov_dstlb_list = [seg2cov.get((gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample), 0.0) for sample in samplesb]

								if any(x >= args.min_dstlCov for x in [np.mean(cov_dstla_list), np.mean(cov_dstlb_list)]) and all (x >= args.min_prxlCov for x in [np.mean(cov_prxla_list), np.mean(cov_prxlb_list)]):

									cov_mean_dstla_scaled, cov_mean_dstlb_scaled, dstl_prxl_a, dstl_prxl_b, cv2a, cv2b, total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered = test_diffl_cp(cov_prxla_list, cov_dstla_list, cov_prxlb_list, cov_dstlb_list,
										total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered)

									if dstl_prxl_a != 'NA' and dstl_prxl_b != 'NA':
										# test_type = 'Tandem_' + test_type
										test_type = get_test_type(strand, side_prxla[0], 'Tandem')
										tested_count += 1
										o5.write('\t'.join(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											gene, test_type, side_prxla,
											str(cstart_prxla) + '-' + str(cend_prxla),
											str(cstart_dstla) + '-' + str(cend_dstla),
											float(ru_prxla) - float(ru_prxlb),
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb,
											cov_mean_dstla_scaled, cov_mean_dstlb_scaled, cv2a, cv2b,
											len(samplesa), len(samplesb)])) + '\n')
								else:
									total_cp_min_dstlCov_filtered += 1
						else:
							test_type = get_test_type(strand, side_prxla[0], 'AE')
							ru_prxlb = 'NA'
							if args.verbose:
								print '-> change point in A, distal or non-consecutive in B' # use distal/proximal coordinates from A only
							if (side_prxla[0] == 'L' and cstart_prxla != cend_dstla) or (side_prxla[0] == 'R' and cend_prxla != cstart_dstla):
								if args.verbose:
									print '-> change point in A, distal in B: not consecutive'
								gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
									test_type, side_prxla,
									str(cstart_prxlb) + '-' + str(cend_prxlb),
									str(cstart_dstlb) + '-' + str(cend_dstlb),
									'NA',
									ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
							else:
								# === test ===
								cov_prxla_list = [seg2cov.get((gene, chrom, str(cstart_prxla), str(cend_prxla), strand, sample), 0.0) for sample in samplesa]
								cov_prxlb_list = [seg2cov.get((gene, chrom, str(cstart_prxla), str(cend_prxla), strand, sample), 0.0) for sample in samplesb]
								cov_dstla_list = [seg2cov.get((gene, chrom, str(cstart_dstla), str(cend_dstla), strand, sample), 0.0) for sample in samplesa]
								cov_dstlb_list = [seg2cov.get((gene, chrom, str(cstart_dstla), str(cend_dstla), strand, sample), 0.0) for sample in samplesb]

								if any(x >= args.min_dstlCov for x in [np.mean(cov_dstla_list), np.mean(cov_dstlb_list)]) and all (x >= args.min_prxlCov for x in [np.mean(cov_prxla_list), np.mean(cov_prxlb_list)]):
									cov_mean_dstla_scaled, cov_mean_dstlb_scaled, dstl_prxl_a, dstl_prxl_b, cv2a, cv2b, total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered = test_diffl_cp(cov_prxla_list, cov_dstla_list, cov_prxlb_list, cov_dstlb_list,
										total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered)

									if dstl_prxl_a != 'NA' and dstl_prxl_b != 'NA':
										# all proximal non-zero & distal < proximal: write output for testing in R
										tested_count += 1
										o5.write('\t'.join(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											gene, test_type, side_prxla,
											str(cstart_prxla) + '-' + str(cend_prxla),
											str(cstart_dstla) + '-' + str(cend_dstla),
											'NA',
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb,
											cov_mean_dstla_scaled, cov_mean_dstlb_scaled, cv2a, cv2b,
											len(samplesa), len(samplesb)])) + '\n')
								else:
									total_cp_min_dstlCov_filtered += 1
					else:
						test_type = get_test_type(strand, side_prxla[0], 'AE')
						ru_prxlb = 'NA'
						ru_dstlb = 'NA'
						if gene in gene2dstlb:
							if (side_prxla[0] == 'L' and cstart_prxla != cend_dstla) or (side_prxla[0] == 'R' and cend_prxla != cstart_dstla):
								if args.verbose:
									print '-> change point in A, not in B: not consectuive ==> distal in A only'
								gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
									test_type, side_prxla,
									str(cstart_prxla) + '-' + str(cend_prxla),
									str(cstart_dstla) + '-' + str(cend_dstla),
									'NA',
									ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
							else:
								if args.verbose:
									print '-> change point in A, not in B: test'
								cov_prxla_list = [seg2cov.get((gene, chrom, str(cstart_prxla), str(cend_prxla), strand, sample), 0.0) for sample in samplesa]
								cov_prxlb_list = [seg2cov.get((gene, chrom, str(cstart_prxla), str(cend_prxla), strand, sample), 0.0) for sample in samplesb]
								cov_dstla_list = [seg2cov.get((gene, chrom, str(cstart_dstla), str(cend_dstla), strand, sample), 0.0) for sample in samplesa]
								cov_dstlb_list = [seg2cov.get((gene, chrom, str(cstart_dstla), str(cend_dstla), strand, sample), 0.0) for sample in samplesb]

								if any(x >= args.min_dstlCov for x in [np.mean(cov_dstla_list), np.mean(cov_dstlb_list)]) and all (x >= args.min_prxlCov for x in [np.mean(cov_prxla_list), np.mean(cov_prxlb_list)]):
									cov_mean_dstla_scaled, cov_mean_dstlb_scaled, dstl_prxl_a, dstl_prxl_b, cv2a, cv2b, total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered = test_diffl_cp(cov_prxla_list, cov_dstla_list, cov_prxlb_list, cov_dstlb_list,
										total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered)

									if dstl_prxl_a != 'NA' and dstl_prxl_b != 'NA':
										tested_count += 1
										o5.write('\t'.join(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											gene, test_type, side_prxla,
											str(cstart_prxla) + '-' + str(cend_prxla),
											str(cstart_dstla) + '-' + str(cend_dstla),
											'NA',
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb,
											cov_mean_dstla_scaled, cov_mean_dstlb_scaled, cv2a, cv2b,
											len(samplesa), len(samplesb)])) + '\n')
									elif dstl_prxl_a != 'NA':
										if args.verbose:
											print '--> print to output separately'
										gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											test_type, side_prxla,
											str(cstart_prxla) + '-' + str(cend_prxla),
											str(cstart_dstla) + '-' + str(cend_dstla),
											'NA',
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
								else:
									total_cp_min_dstlCov_filtered += 1
						elif args.verbose:
							print '-> change point in A, not in B: gene not expressed or all change points were noise in B. DO NOT OUTPUT'
				else:
					ru_prxla = 'NA'
					if (cp_start, cp_end, label, gene) in cp2dstlb:
						if (cp_start, cp_end, label, gene) in cp2prxlb:
							test_type = get_test_type(strand, side_prxlb[0], 'AE')
							if (side_prxlb[0] == 'L' and cstart_prxlb != cend_dstlb) or (side_prxlb[0] == 'R' and cend_prxlb != cstart_dstlb):
								if args.verbose:
									print '-> change point in B, distal in A: not consecutive'
								gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
									test_type, side_prxlb,
									str(cstart_prxlb) + '-' + str(cend_prxlb),
									str(cstart_dstlb) + '-' + str(cend_dstlb),
									'NA',
									ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
							else:
								if args.verbose:
									print '-> change point in B, distal in A' # use distal/proximal coordinates from B only
								# === test ===
								cov_prxla_list = [seg2cov.get((gene, chrom, str(cstart_prxlb), str(cend_prxlb), strand, sample), 0.0) for sample in samplesa]
								cov_prxlb_list = [seg2cov.get((gene, chrom, str(cstart_prxlb), str(cend_prxlb), strand, sample), 0.0) for sample in samplesb]
								cov_dstla_list = [seg2cov.get((gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample), 0.0) for sample in samplesa]
								cov_dstlb_list = [seg2cov.get((gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample), 0.0) for sample in samplesb]

								for sample in samplesa:
									if (gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample) in seg2cov:
										print '-> cov a:', gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample, seg2cov[(gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample)]
									else:
										print '-> none a!:', gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample

								for sample in samplesb:
									if (gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample) in seg2cov:
										print '-> cov b:', gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample, seg2cov[(gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample)]
									else:
										print '-> none b!:', gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample

								if any(x >= args.min_dstlCov for x in [np.mean(cov_dstla_list), np.mean(cov_dstlb_list)]) and all (x >= args.min_prxlCov for x in [np.mean(cov_prxla_list), np.mean(cov_prxlb_list)]):
									cov_mean_dstla_scaled, cov_mean_dstlb_scaled, dstl_prxl_a, dstl_prxl_b, cv2a, cv2b, total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered = test_diffl_cp(cov_prxla_list, cov_dstla_list, cov_prxlb_list, cov_dstlb_list,
										total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered)

									if dstl_prxl_a != 'NA' and dstl_prxl_b != 'NA':
										tested_count += 1
										o5.write('\t'.join(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											gene, test_type, side_prxlb,
											str(cstart_prxlb) + '-' + str(cend_prxlb),
											str(cstart_dstlb) + '-' + str(cend_dstlb),
											'NA',
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb,
											cov_mean_dstla_scaled, cov_mean_dstlb_scaled, cv2a, cv2b,
											len(samplesa), len(samplesb)])) + '\n')
								else:
									total_cp_min_dstlCov_filtered += 1
						elif args.verbose:
							print '-> distal in both A & B: not differential: DO NOT OUTPUT (there may be other tandem change points proximal of this distal point'
					else:
						ru_prxlb = 'NA'
						ru_dstlb = 'NA'
						test_type = get_test_type(strand, side_dstla[0], 'AE')
						if gene in gene2dstlb:
							if args.verbose:
								print '--> distal only in A, not in B: change point distal in one condition but not the other: output'
							dstl_prxl_a = 'NA'
							dstl_prxl_b = 'NA'
							gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
								test_type, side_dstla,
								'NA',
								str(cstart_dstla) + '-' + str(cend_dstla),
								'NA',
								ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
						elif args.verbose:
							print '--> distal in A, not in B: gene not expressed or all change points were noise in B. DO NOT OUTPUT'
			else:
				ru_prxla = 'NA'
				ru_dstla = 'NA'
				if (cp_start, cp_end, label, gene) in cp2dstlb:
					if (cp_start, cp_end, label, gene) in cp2prxlb:
						test_type = get_test_type(strand, side_prxlb[0], 'AE')
						if gene in gene2dstla:
							if (side_prxlb[0] == 'L' and cstart_prxlb != cend_dstlb) or (side_prxlb[0] == 'R' and cend_prxlb != cstart_dstlb):
								if args.verbose:
									print '-> not in A, change point in B: not consecutive ==> distal in B only'
								dstl_prxl_b = 'NA'
								dstl_prxl_a = 'NA'
								gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
									test_type, side_prxlb,
									str(cstart_prxlb) + '-' + str(cend_prxlb),
									str(cstart_dstlb) + '-' + str(cend_dstlb),
									'NA',
									ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
							else:
								if args.verbose:
									print '-> not in A, change point in B: test'
								# === test ===
								cov_prxla_list = [seg2cov.get((gene, chrom, str(cstart_prxlb), str(cend_prxlb), strand, sample), 0.0) for sample in samplesa]
								cov_prxlb_list = [seg2cov.get((gene, chrom, str(cstart_prxlb), str(cend_prxlb), strand, sample), 0.0) for sample in samplesb]
								cov_dstla_list = [seg2cov.get((gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample), 0.0) for sample in samplesa]
								cov_dstlb_list = [seg2cov.get((gene, chrom, str(cstart_dstlb), str(cend_dstlb), strand, sample), 0.0) for sample in samplesb]

								if any(x >= args.min_dstlCov for x in [np.mean(cov_dstla_list), np.mean(cov_dstlb_list)]) and all (x >= args.min_prxlCov for x in [np.mean(cov_prxla_list), np.mean(cov_prxlb_list)]):
									cov_mean_dstla_scaled, cov_mean_dstlb_scaled, dstl_prxl_a, dstl_prxl_b, cv2a, cv2b, total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered = test_diffl_cp(cov_prxla_list, cov_dstla_list, cov_prxlb_list, cov_dstlb_list,
										total_cp_min_dstlCov, total_cp_nonzero_prxl, total_cp_nonzero_prxl_filtered, total_cp_strictly_decreasing, total_cp_strictly_decreasing_filtered)

									if dstl_prxl_a != 'NA' and dstl_prxl_b != 'NA':
										tested_count += 1
										o5.write('\t'.join(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											gene, test_type, side_prxlb,
											str(cstart_prxlb) + '-' + str(cend_prxlb),
											str(cstart_dstlb) + '-' + str(cend_dstlb),
											'NA',
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb,
											cov_mean_dstla_scaled, cov_mean_dstlb_scaled, cv2a, cv2b,
											len(samplesa), len(samplesb)])) + '\n')

									elif dstl_prxl_b != 'NA':
										gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
											test_type, side_prxlb,
											str(cstart_prxlb) + '-' + str(cend_prxlb),
											str(cstart_dstlb) + '-' + str(cend_dstlb),
											'NA',
											ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
						elif args.verbose:
							print '--> not in A, change point in B: gene not expressed or all change points were noise in A. DO NOT OUTPUT'
					else:
						ru_prxlb = 'NA'
						if gene in gene2dstla:
							if args.verbose:
								print '--> not in A, distal only in B: change point distal in one condition but not the other: output (e.g. HEPH)'
							dstl_prxl_b = 'NA'
							dstl_prxl_a = 'NA'
							test_type = get_test_type(strand, side_dstlb[0], 'AE')
							gene2diffcp[gene].append(map(str, [':'.join(map(str, [label, chrom, cp_start, cp_end])),
								test_type, side_dstlb,
								'NA',
								str(cstart_dstlb) + '-' + str(cend_dstlb),
								'NA',
								ru_prxla, ru_prxlb, ru_dstla, ru_dstlb]))
						elif args.verbose:
							print '--> not in A, distal only in B: gene not expressed or all change points were noise in B. DO NOT OUTPUT'
				else:
					print '-> not in A or B'
					sys.exit(1)
		o5.close()

		# run differential expression test
		if tested_count == 0:
			print '-> no change points tested!'
		else:
			print '- differential expression test:'
			rscript = os.path.join(os.path.dirname(functions.__file__), 'diff_test.R')
			cmd = map(str, ['Rscript', '--vanilla', rscript, file_test, args.output, args.pmax, args.dtop_abs_dif_min])
			stdout, stderr = run_command(cmd)
			print stderr
			print stdout

			total_pmax = ''
			total_pmax_difmin = ''
			test_lines = stdout.split('\n')
			for test_line in test_lines:
				if 'p <= ' + str(args.pmax) + '"' in test_line:
					total_pmax = test_line.split(' ')[1].replace('"', '')
				if 'p <= ' + str(args.pmax) + ' & abs(RUa - RUb) >= ' + str(args.dtop_abs_dif_min) + '"' in test_line:
					total_pmax_difmin = test_line.split(' ')[1].replace('"', '')
			if total_pmax == '' or total_pmax_difmin == '':
				sys.stderr.write('EXIT: diff_test.R did not execute properly\n')
				sys.exit(1)

			# === write output: ===
			# totals
			o8 = open(file_totals, 'w')
			o8.write('\t'.join(['total cp',
				'total cp with distal and proximal segment reads/bp avg across samples >= ' +
				str(args.min_dstlCov) + ' and >= ' + str(args.min_prxlCov) + ' in >= 1 condition respectively',
				'total cp strictly decreasing from proximal to distal',
				'p <= ' + str(args.pmax), 'p <= ' + str(args.pmax) + ' & abs(RUa - RUb) >= ' + str(args.dtop_abs_dif_min)]) + '\n')
			o8.write('\t'.join(map(str, [total_cp, total_cp_min_dstlCov,
				total_cp_strictly_decreasing,
				total_pmax, total_pmax_difmin])) + '\n')

			# bed files
			o6 = open(outfile_cp_allTested, 'w')
			o6.write('track name=\"' + os.path.basename(args.output) + ' all tested\"\n')
			o7 = open(outfile_cp_diff, 'w')
			o7.write('track name=\"' + os.path.basename(args.output) + ' differential\"\n')

			# tandem change points
			dtop_abs_dif_min_count = 0
			landmarker_test_file = args.output + '.txt'
			with open(landmarker_test_file, 'r') as f:
				for line in f:
					if not line.startswith('cp'):
						(cp, gene, test_type, side, seg_prxl, seg_dstl, ru_dif, ru_prxla, ru_prxlb, ru_dstla, ru_dstlb, pval, pBH) = line.rstrip().split('\t')
						(cplabel, chrom, cpstart, cpend) = cp.split(':')
						(gs, gchrom, gstart, gend, strand) = gene.split(':')
						o6.write('\t'.join([chrom, cpstart, cpend,
							';'.join([cplabel, test_type, side,
								':'.join([gs, gstart, gend]),
								seg_prxl, seg_dstl, ru_dif,
								ru_prxla, ru_prxlb,
								pBH]),
							'0', strand]) + '\n')
						if ru_dif != 'NA' and pBH != 'NA':
							if abs(float(ru_dif)) >= args.dtop_abs_dif_min and float(pBH) <= args.pmax:
								dtop_abs_dif_min_count += 1
								o7.write('\t'.join([chrom, cpstart, cpend,
									';'.join([cplabel, test_type, side,
										':'.join([gs, gstart, gend]),
										seg_prxl, seg_dstl, ru_dif,
										ru_prxla, ru_prxlb,
										pBH]),
									'0', strand]) + '\n')

			# condition-specific change points
			for gene in gene2diffcp:
				cp_list = gene2diffcp[gene]
				for each_cp in cp_list:
					(cp, test_type, side, seg_prxl, seg_dstl, ru_dif, ru_prxla, ru_prxlb, ru_dstla, ru_dstlb) = each_cp
					(cplabel, chrom, cpstart, cpend) = cp.split(':')
					(gs, gchrom, gstart, gend, strand) = gene.split(':')

					o6.write('\t'.join([chrom, cpstart, cpend,
						';'.join([cplabel, test_type, side,
							':'.join([gs, gstart, gend]),
							seg_prxl, seg_dstl, ru_dif, ru_prxla, ru_prxlb,
							'NA']),
						'0', strand]) + '\n')

					if ru_dif != 'NA':
						if abs(float(ru_dif)) >= args.dtop_abs_dif_min:
							dtop_abs_dif_min_count += 1
							o7.write('\t'.join([chrom, cpstart, cpend,
								';'.join([cplabel, test_type, side,
									':'.join([gs, gstart, gend]),
									seg_prxl, seg_dstl, ru_dif, ru_prxla, ru_prxlb,
									'NA']),
								'0', strand]) + '\n')
			o6.close()
			o7.close()

			print dtop_abs_dif_min_count, 'total cp with abs(RUa - RUb) >=', args.dtop_abs_dif_min
	else:
		sys.stderr.write('EXIT: please input two conditions\n')
		sys.exit(1)

	# remove intermediate files
	if not args.keep:
		for file in [landmarker_test_file, file_test]:
			os.remove(file)

	print 'finished:', str(datetime.now().time())

# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
