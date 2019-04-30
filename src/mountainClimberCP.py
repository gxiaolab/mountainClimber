#!/user/bin/python -tt
"""
Identify change points throughout the entire TU of each individual sample.
The premise of this approach is to identify significant change points in the
cumulative read sum (CRS) distribution as a function of position. It identifies
the following change point types: DistalTSS, TandemTSS, DistalPolyA, TandemAPA,
Junction, Exon, and Intron.
"""


import os
import sys
import argparse
from datetime import datetime
from scipy import stats 			# v0.15.1
import math
import numpy as np 					# v1.10.4
import peakutils					# v1.0.3
from collections import defaultdict, OrderedDict, deque
import pysam 						# v0.9.0
import pybedtools as pb
from functions import sort_bedfile

# median filter
from bisect import bisect_left, insort
from itertools import islice

# plotting
import matplotlib
matplotlib.use('pdf')  # force matplotlib to not use any Xwindows backend for plotting
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages


def bedgraph_per_gene_ss(genes, bg_plus, bg_minus, bgfile):
	"""
	bedtools intersect genes with each of bg_plus and bg_minus.
	Run separately so that gene coverage is consecutive by strand.
	"""
	# === split annotation ===
	plus_bed = bgfile + '.genes.plus'
	minus_bed = bgfile + '.genes.minus'
	p = open(plus_bed, 'w')
	m = open(minus_bed, 'w')
	with open(genes, 'r') as f:
		for line in f:
			if not line.startswith('track'):
				strand = line.rstrip().split('\t')[5]
				if strand == '+':
					p.write(line)
				elif strand == '-':
					m.write(line)
				else:
					sys.stderr.write('do not recognize strand: ' + strand + '\n')
					sys.stderr.write(line)
					sys.exit(1)
	p.close()
	m.close()

	# === bedtools intersect: concatenate + & - strands ===
	pb.BedTool(plus_bed).intersect(bg_plus, wo=True, sorted=True).saveas(bgfile + '.plus')
	pb.BedTool(minus_bed).intersect(bg_minus, wo=True, sorted=True).saveas(bgfile + '.minus')

	t = open(bgfile, 'w')
	t.write(open(bgfile + '.plus').read())
	t.write(open(bgfile + '.minus').read())
	t.close()

	for file in [bgfile + '.plus', bgfile + '.minus', plus_bed, minus_bed]:
		os.remove(file)


def bedgraph_per_gene_nss(genes, bg, bgfile):
	"""Bedtools intersect, non-strand-specific"""
	pb.BedTool(genes).intersect(bg, wo=True, sorted=True).saveas(bgfile)


def get_exon_cov(exon_list, cov_list):
	"""Calculate average reads/bp in each exon and report the maximum exon coverage."""
	total_sum = 0
	total_len = 0
	max_cov = 0
	for e, exon in enumerate(exon_list):
		this_start, this_end = map(int, exon.split(':'))
		this_cov = cov_list[this_start:this_end]
		total_sum += sum(this_cov)
		total_len += this_end - this_start
		this_exon_cov = float(sum(this_cov)) / float(this_end - this_start)
		if this_exon_cov > max_cov:
			max_cov = this_exon_cov

	cov_avg_all_exons = float(total_sum) / float(total_len)
	return cov_avg_all_exons, max_cov


def crs(cov_array):
	"""Calculate cumulative read sum """
	vert_array = np.insert(np.ediff1d(cov_array), [0], 0)
	vert_sum_array = np.cumsum(np.absolute(vert_array))
	if max(vert_sum_array) == 0:
		vert_sum_norm_array = ['NA']
	else:
		vert_sum_norm_array = vert_sum_array / max(vert_sum_array)
	return vert_sum_norm_array, vert_array


def ks_test(vert_sum_array, make_plots, out_prefix):
	"""KS test: cumulative distance vs. line y=ax"""
	line_array = np.arange(0, max(vert_sum_array), max(vert_sum_array) / vert_sum_array.size)
	ks_stat, ksp = stats.ks_2samp(vert_sum_array, line_array)

	y0 = vert_sum_array[0]
	xmax = vert_sum_array.size - 1
	ymax = max(vert_sum_array)
	slope = (ymax - y0) / xmax

	if slope == 0:
		ksp = 1

	if make_plots:
		x1 = np.linspace(0, xmax, vert_sum_array.size) / 1000
		y1 = vert_sum_array.tolist()
		x2 = [0, float(xmax) / float(1000)]
		y2 = [y0, ymax]

		out_plot = out_prefix + '_crs.pdf'
		pdf = PdfPages(out_plot)
		fig = pyplot.figure(figsize=(3, 3))
		pyplot.plot(x1, y1, color='k')
		pyplot.plot(x2, y2, color='0.75')
		pyplot.title('KS test p = ' + str(round(ksp, 3)))
		pyplot.xlabel('Position (kb)')
		pyplot.ylabel('Cumulative Vertical Distance')
		pyplot.gcf().subplots_adjust(bottom=0.3, left=0.3)
		pdf.savefig()
		pdf.close()
		pyplot.close(fig)

	return ksp


def distance(x1, y1, x2, y2, x0, y0):
	"""Calculate distance from a point x0, y0 to a line defined by x1, y1 and x2, y2"""
	num = (y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1
	den = math.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)
	return float(num) / float(den)


def get_slope(x1, y1, x2, y2):
	""""Get slope between 2 points"""
	num = y2 - y1
	den = x2 - x1
	return float(num) / float(den)


def median_filter(seq, window_size):
	"""
	Median filter.
	Reference: https://groups.google.com/forum/#!topic/comp.lang.python/0OARyHF0wtA
	"""
	zeros = [0] * (window_size / 2)  # pad zeros
	seq = zeros + seq + zeros

	seq = iter(seq)
	d = deque()
	s = []
	result = []

	# first window
	for item in islice(seq, window_size):
		d.append(item)
		insort(s, item)
	if len(s) % 2 == 0:
		med = float((s[len(d) / 2 - 1] + s[len(d) / 2])) / float(2)  # even number: take average of 2
		result.append(med)
	else:
		result.append(float(s[len(d) / 2]))

	# all other windows
	m = window_size / 2
	for item in seq:
		old = d.popleft()
		d.append(item)
		del s[bisect_left(s, old)]
		insort(s, item)
		if len(s) % 2 == 0:
			med = float((s[m - 1] + s[m])) / float(2)  # even number: take average of 2
			result.append(med)
		else:
			result.append(float(s[m]))
	return result


def get_linedist(vert_sum_array, denoise_winsize):
	"""De-noise with median filter & get distance of all points to the line y=ax"""
	x1 = 0
	y1 = 0
	x2 = vert_sum_array.size - 1
	y2 = 1
	x0 = np.arange(vert_sum_array.size)
	y0 = vert_sum_array

	# distance from (x0, y0) to the line formed by (x1, y1) and (x2, y2)
	line_dist_array = ((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / math.sqrt(
		(y2 - y1)**2 + (x2 - x1)**2)

	# de-noising: median filter. window size should be odd for medfilt to function
	denoise_winsize = denoise_winsize + 1 if denoise_winsize % 2 == 0 else denoise_winsize
	line_dist_list_denoise = median_filter(line_dist_array.tolist(), denoise_winsize)
	line_dist_array_denoise = np.asarray(line_dist_list_denoise)

	return line_dist_array, line_dist_array_denoise


def get_introns_exons(intron2jxnCount, gene_len, chrom, start, end, strand):
	"""get introns (0-based): junction read-spanning regions with at least 2 reads"""
	jxn_list = []
	start2end = {}
	end2start = {}
	either2either = {}
	(start, end) = map(int, [start, end])
	for intron in intron2jxnCount:
		x = intron.split(':')
		if len(x) == 4:
			[ichrom, istart, iend, istrand] = x
		else:
			[ichrom, istart, iend] = x
			istrand = 0
		if ichrom == chrom and istrand == strand:
			# start is 0-based, last pos on exon; end is 1-based, first pos on following exon
			[istart, iend] = map(int, [istart, iend])
			if int(start) <= istart < int(end) and int(start) < iend <= int(end):  # intron in this gene
				if istart - start != 0:
					jxn_list.append(istart - start)
				if iend - 1 - start != gene_len:
					jxn_list.append(iend - 1 - start)  # 1-based -> 0-based

				# get exons
				if istart not in start2end:
					start2end[istart] = iend
				elif iend - istart < start2end[istart] - istart:
					start2end[istart] = iend
				if iend not in end2start:
					end2start[iend] = istart
				elif iend - istart < iend - end2start[iend]:
					end2start[iend] = istart

	# get exon list
	intron_list = []
	for istart in start2end:
		intron_list.append((istart, start2end[istart]))
	for iend in end2start:
		if (end2start[iend], iend) not in intron_list:
			intron_list.append((end2start[iend], iend))

	intron_list = sorted(intron_list)
	exon_list = []
	for i, intron in enumerate(intron_list):
		estart = intron_list[i - 1][1] - 1 - start if i > 0 else 0
		eend = intron[0] - start
		if eend > estart:
			exon_list.append(':'.join(map(str, [estart, eend])))
		if i == len(intron_list) - 1:  # last exon
			estart = intron[1] - 1 - start
			eend = gene_len
			exon_list.append(':'.join(map(str, [estart, eend])))
	if len(exon_list) == 0:
		exon_list.append(':'.join(map(str, [0, end - 1 - start])))

	return sorted(list(set(jxn_list))), exon_list, intron_list


def increase_precision(peak_inds, denoise_winsize, line_dist_array, max_or_min, peak_thresh):
	"""call peak in pre-de-noised data"""
	for i, ind in enumerate(peak_inds):
		start = max(ind - denoise_winsize, 0)
		end = min(ind + denoise_winsize, len(line_dist_array))
		window = line_dist_array[start:end]

		if max_or_min == 'max':
			peak_inds_new = peakutils.peak.indexes(window, thres=peak_thresh,
												   min_dist=denoise_winsize * 2)
		elif max_or_min == 'min':
			peak_inds_new = peakutils.peak.indexes(-1 * window, thres=peak_thresh,
												   min_dist=denoise_winsize * 2)
		else:
			sys.exit(1)

		if len(peak_inds_new) == 1:
			peak_ind_new = peak_inds_new[0] + start
			peak_inds[i] = peak_ind_new
		elif len(peak_inds_new) != 0:
			sys.stderr.write('too many new peak indexes detected')
			sys.stderr.write(' '.join([start, end, denoise_winsize, num_mins]) + '\n')
			sys.exit(1)

	return np.sort(peak_inds)


def get_ttest(peak_inds, cov_array, denoise_winsize, test_thresh):
	"""Calculate t-test of read coverage per nucleotide before vs. after the change point"""
	peak_inds_ttest = []
	ind2tp = {}
	if len(peak_inds) > 0:
		for ind in peak_inds:
			next_cut = min(cov_array.size - 1, ind + denoise_winsize)
			prev_cut = max(0, ind - denoise_winsize)
			if ind - prev_cut > 1 and next_cut - ind > 1:
				cov_before = cov_array[prev_cut:ind]
				cov_after = cov_array[ind:next_cut]
				if np.std(cov_before) != 0.0 or np.std(cov_after) != 0.0:  # make sure denominator != 0
					(t_stat, tp) = stats.ttest_ind(cov_before, cov_after, equal_var=False)
				else:
					tp = 1

				if tp < test_thresh:
					peak_inds_ttest.append(ind)
					ind2tp[ind] = tp

	peak_inds_ttest = np.sort(peak_inds_ttest).tolist()
	return peak_inds_ttest, ind2tp


def get_end_cov(utr5, utr3, cov_list, peak_inds_ttest):
	"""Require increasing & decreasing coverage at left & right ends respecitvely"""
	to_delete = []
	if len(utr5) > 0:
		u2prev_mean = {}
		for u, u5cp in enumerate(utr5):
			if u < len(utr5) - 1:
				prev_cp = 0 if u == 0 else prev_cp
				next_cp = utr5[u + 1]
				this_mean = np.mean(cov_list[u5cp:next_cp])
				prev_mean = np.mean(cov_list[prev_cp:u5cp])
				u2prev_mean[u] = prev_mean
				del_ind = u
				prev_u = u - 1
				flag = 0
				while this_mean < prev_mean:
					flag = 1
					to_delete.append(del_ind)
					if prev_u > 0:
						if utr5[prev_u] not in to_delete:
							this_mean = np.mean(cov_list[utr5[prev_u]:next_cp])
							prev_mean = u2prev_mean[prev_u]
							del_ind = prev_u
							prev_cp = utr5[prev_u]
							prev_u -= 1
						else:
							prev_u -= 1
					else:
						break
				if flag != 1:
					prev_cp = u5cp
			else:
				prev_cp = u5cp

	if len(utr3) > 0:
		u2prev_mean = {}
		for u, u3cp in enumerate(utr3):
			if u > 0:
				next_cp = len(cov_list) if u == len(utr3) - 1 else utr3[u + 1]
				this_mean = np.mean(cov_list[u3cp:next_cp])
				prev_mean = np.mean(cov_list[prev_cp:u3cp])
				u2prev_mean[u] = prev_mean
				del_ind = len(peak_inds_ttest) - len(utr3) + u
				prev_u = u - 1
				flag = 0
				while this_mean > prev_mean:
					flag = 1
					to_delete.append(del_ind)
					if prev_u > 0:
						if utr3[prev_u] not in to_delete:
							this_mean = np.mean(cov_list[utr3[prev_u]:next_cp])
							prev_mean = u2prev_mean[prev_u]
							del_ind = len(peak_inds_ttest) - len(utr3) + prev_u
							prev_cp = utr3[prev_u]
							prev_u -= 1
						else:
							prev_u -= 1
					else:
						break
				if flag != 1:
					prev_cp = u3cp
			else:
				prev_cp = u3cp
	peak_inds_ttest = sorted(np.delete(peak_inds_ttest, to_delete))
	return peak_inds_ttest


def infer_strand(intron_list, chrom, genome, verbose):
	"""Infer strand of TU from non-strand-specific RNA-Seq by checking the GT-AG splice site signal"""
	count_minus = 0
	count_plus = 0
	count_unknown = 0

	intron_label_list = []
	for (istart, iend) in intron_list:
		# beginning of intron (junction is 0-based & 1st position on intron)
		start_coord = chrom + ':' + str(istart + 1) + '-' + str(istart + 2)
		start_seq = pysam.faidx(genome, start_coord)

		# end of intron (junction is 1-based & last position on intron)
		end_coord = chrom + ':' + str(iend - 1) + '-' + str(iend)
		end_seq = pysam.faidx(genome, end_coord)

		if start_seq.split('\n')[1].upper() == 'CT': 	# 3' splice site AG on minus strand
			count_minus += 1
			intron_label_list.append('-')
		elif start_seq.split('\n')[1].upper() == 'GT': 	# 5' splice site GT on plus strand
			count_plus += 1
			intron_label_list.append('+')
		else:
			count_unknown += 1
			intron_label_list.append('?')
		if end_seq.split('\n')[1].upper() == 'AC':		# 5' splice site GT on minus strand
			count_minus += 1
			intron_label_list.append('-')
		elif end_seq.split('\n')[1].upper() == 'AG':		# 3' splice site AG on plus strand
			count_plus += 1
			intron_label_list.append('+')
		else:
			count_unknown += 1
			intron_label_list.append('?')

	if max(count_minus, count_plus, count_unknown) == count_unknown:
		strand_inferred = 'NA'
	elif count_minus == count_plus:
		strand_inferred = 'NA'
	elif max(count_minus, count_plus, count_unknown) == count_minus:
		strand_inferred = '-'
	elif max(count_minus, count_plus, count_unknown) == count_plus:
		strand_inferred = '+'

	if verbose:
		print '- inferred strand:', strand_inferred
		print '  - counts: minus', count_minus, 'plus', count_plus, 'unknown', count_unknown

	return strand_inferred


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									 description='Identify change points throughout each TU using the cumulative read sum (CRS).')
	group = parser.add_argument_group('Input')
	group.add_argument('-i', '--input_bg', dest='input_bg', type=str, metavar='',
					   help='Bedgraph: non-strand-specific or plus strand.')
	group.add_argument('-m', '--input_bg_minus', dest='input_bg_minus', type=str, metavar='',
					   help='Bedgraph, minus strand (strand-specific only).')
	group.add_argument('-g', '--input_regions', dest='input_regions', type=str, metavar='',
					   help='Bed file of transcription units.')
	group.add_argument('-j', '--junc', dest='junc', type=str, metavar='',
					   help='Bed file of junction read counts.')
	group.add_argument('-x', '--genome', dest='genome', type=str, metavar='',
					   help='Genome fasta file.')

	group = parser.add_argument_group('Parameters')
	group.add_argument('-a', '--peak_thresh', dest='peak_thresh', type=float, default=-1.0, metavar='',
					   help='Normalized threshold (float between [0., 1.]). Only the peaks with amplitude higher than the threshold will be detected. -1.0 indicates optimizing between [0.01, 0.05, 0.1] for each TU.')
	group.add_argument('-d', '--peak_min_dist', dest='peak_min_dist', type=int, default=-1, metavar='',
					   help='Minimum distance between each detected peak. The peak with the highest amplitude is preferred to satisfy this constraint. -1 indicates optimizing between [10, 50] for each TU.')
	group.add_argument('-w', '--winsize', dest='winsize', type=int, default=-1, metavar='',
					   help='Window size for de-noising and increasing precision. -1 indicates optimizing between [50, max(100, gene_length / 100) * 2].')
	group.add_argument('-t', '--test_thresh', dest='test_thresh', type=float, default=0.001, metavar='',
					   help='Maximum p-value threshold for KS test and t-test.')
	group.add_argument('-l', '--min_length', dest='min_length', type=int, default=1000, metavar='',
					   help='Minimum gene length for running mountain climber.')
	group.add_argument('-e', '--min_expn', dest='min_expn', type=int, default=10, metavar='',
					   help='Minimum expression level (average # reads per bp)')
	group.add_argument('-s', '--min_expn_distal', dest='min_expn_distal', type=int, default=1, metavar='',
					   help='Minimum distal expression level (average # reads per bp).')
	group.add_argument('-f', '--fcthresh', dest='fcthresh', type=float, default=1.5, metavar='',
					   help='Minimum fold change.')
	group.add_argument('-u', '--juncdist', dest='juncdist', type=int, default=10, metavar='',
					   help='Minimum distance to exon-intron junction.')
	group.add_argument('-n', '--minjxncount', dest='minjxncount', type=int, default=2, metavar='',
					   help='Minimum junction read count.')
	group.add_argument('-z', '--max_end_ru', dest='max_end_ru', type=float, default=0.01, metavar='',
					   help='Maximum end relative usage = coverage of end / max segment coverage.')

	group = parser.add_argument_group('Output')
	group.add_argument('-o', '--output', dest='output', type=str, metavar='',
					   help='Output prefix. Bed file of change points has name field = CPlabel:gene:TUstart:TUend:inferred_strand:winsize:peak_thresh:peak_min_dist:segment1coverage:segment2coverage:segment1sd:segment2sd:exon_coverage. Score = log2(fold change).')
	group.add_argument('-p', '--plot', dest='plot', action='store_true',
					   help='Plot the cumulative read sum (CRS), the distance from CRS to line y=ax, and the coverage with predicted change points.')
	group.add_argument('-v', '--verbose', dest='verbose', action='store_true',
					   help='Print progress.')
	args = parser.parse_args()

	# --------------------------------------------------
	# main routine
	# --------------------------------------------------
	print '\njob starting:', str(datetime.now().time())

	if not args.input_bg:
		sys.stderr.write('EXIT: Please provide --input_bg\n')
		sys.exit(1)

	if not args.input_regions:
		sys.stderr.write('EXIT: Please provide --input_regions\n')
		sys.exit(1)

	if not args.output:
		sys.stderr.write('EXIT: Please provide --output\n')
		sys.exit(1)

	if not args.input_bg_minus and not args.genome:
		sys.stderr.write('EXIT: for non-strand-specific RNA-Seq, please provide a --genome file\n')
		sys.exit(1)

	plot_dir = os.path.splitext(args.output)[0] + '_plots'
	if args.plot:
		if os.path.exists(plot_dir):
			for filename in os.listdir(plot_dir):
				os.remove(os.path.join(plot_dir, filename))
		else:
			os.mkdir(plot_dir)

	# === get bedgraph per gene ===
	print '- bedtools intersect', str(datetime.now().time())
	bgfile = os.path.splitext(args.output)[0] + '_gene_bedgraph_intersect.txt'
	if args.input_bg_minus: 	# strand-specific
		bedgraph_per_gene_ss(args.input_regions, args.input_bg, args.input_bg_minus, bgfile)
	else: 						# non-strand-specific
		bedgraph_per_gene_nss(args.input_regions, args.input_bg, bgfile)

	# === get introns from junction counts ===
	print '- getting junction reads', str(datetime.now().time())
	intron2jxnCount = OrderedDict()
	f = open(args.junc, 'r')
	for line in f:
		if not line.startswith('track'):
			x = line.rstrip().split('\t')
			if len(x) == 6:
				(chrom, start, end, name, count, strand) = x
				intron = ':'.join([chrom, start, end, strand])
			elif len(x) == 5:
				(chrom, start, end, name, count) = x
				intron = ':'.join([chrom, start, end])
			else:
				sys.stderr.write('EXIT: did not recognize junction file format\n')
				sys.exit(1)

			if int(count) >= args.minjxncount:
				if intron not in intron2jxnCount:
					intron2jxnCount[intron] = count
				else:
					sys.stderr.write('seen intron')
					sys.exit(1)
	f.close()

	# === find change points in each gene ===
	count_genes_with_reads = 0
	count_min_length_filter = 0
	count_min_length_keep = 0
	count_min_expn_filter = 0
	count_min_expn_keep = 0
	count_high_ksp = 0
	count_vert_sum_max0 = 0
	count_vert2_sum_max0 = 0
	count_no_cps_afterfiltering = 0
	count_no_cps_ttest0 = 0
	count_cps_called = 0
	count_filter123 = 0

	count_genes_with_reads_annotated = 0
	count_min_length_filter_annotated = 0
	count_min_expn_filter_annotated = 0

	print '- finding change points in each gene', str(datetime.now().time())
	maxl = 0
	with open(bgfile, 'r') as f:
		for l, line in enumerate(f):
			maxl = l

	o = open(args.output + '.temp', 'w')
	with open(bgfile, 'r') as f:
		for l, line in enumerate(f):
			# not EOF -> read the line
			if line != '':
				x = line.rstrip().split('\t')

				if len(x) == 11:
					(achrom, astart, aend, ageneid, ascore, astrand, bchrom, bstart, bend, bcov, overlap_len) = x
				elif len(x) == 9:
					(achrom, astart, aend, ageneid, bchrom, bstart, bend, bcov, overlap_len) = x
				elif len(x) == 10:
					(achrom, astart, aend, ageneid, ascore, bchrom, bstart, bend, bcov, overlap_len) = x
				else:
					sys.stderr.write('EXIT: do not recognize bedgraph intersect format\n')
					sys.stderr.write(line)
					sys.exit(1)

				if not args.input_bg_minus:
					astrand = 0

				astart = int(astart)
				aend = int(aend)
				bstart = int(bstart)
				bend = int(bend)
				bcov = float(bcov)
			else:
				x = ''

			if l == 0:  # first line
				prev_gene = ':'.join(x[:5]) if astrand == 0 else ':'.join(x[:6])
				this_start = max(astart, bstart)
				this_end = min(aend, bend)
				prev_cov_array = np.zeros(aend - astart)
				prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

				# === next round ===
				prev_astart = astart
				prev_aend = aend
				prev_geneid = ':'.join(map(str, [ageneid, astart, aend, achrom, astrand]))
				prev_chrom = achrom
				prev_strand = astrand

				if l == maxl:
					cov_sum = sum(prev_cov_array)
					new_start = np.nonzero(prev_cov_array)[0][0]
					new_end = np.nonzero(prev_cov_array)[0][-1]
					prev_cov_array = prev_cov_array[new_start:(new_end + 1)]

					geneid = prev_geneid
					start = prev_astart
					end = prev_aend
					chrom = prev_chrom
					strand = prev_strand
			else:
				this_gene = ':'.join(x[:5]) if astrand == 0 else ':'.join(x[:6])
				if line == '' and this_gene == prev_gene and this_gene == '':  # EOF
					break
				elif this_gene == prev_gene and l != maxl:  # get coverage
					this_start = max(astart, bstart)
					this_end = min(aend, bend)
					prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

					# === next round ===
					prev_astart = astart
					prev_aend = aend
					prev_geneid = ':'.join(map(str, [ageneid, astart, aend, achrom, astrand]))
					prev_chrom = achrom
					prev_strand = astrand
				else:  # finished reading all info for one gene -> call change points
					# === get per-base coverage ===
					if l == maxl:
						# === last line of this gene & last line of the file ===
						this_start = max(astart, bstart)
						this_end = min(aend, bend)
						prev_cov_array[(this_start - astart):(this_end - astart)] += bcov
						cov_sum = sum(prev_cov_array)
						new_start = np.nonzero(prev_cov_array)[0][0]
						new_end = np.nonzero(prev_cov_array)[0][-1]
						prev_cov_array = prev_cov_array[new_start:(new_end + 1)]

						# === call change points in this gene ===
						geneid = prev_geneid
						start = prev_astart
						end = prev_aend
						chrom = prev_chrom
						strand = prev_strand
						cov_array = prev_cov_array
						new_start = new_start + start
						new_end = new_end + start + 1
					else:
						cov_sum = sum(prev_cov_array)
						new_start = np.nonzero(prev_cov_array)[0][0]
						new_end = np.nonzero(prev_cov_array)[0][-1]
						prev_cov_array = prev_cov_array[new_start:(new_end + 1)]

						# === call change points in the previous gene ===
						geneid = prev_geneid
						start = prev_astart
						end = prev_aend
						chrom = prev_chrom
						strand = prev_strand
						cov_array = prev_cov_array
						new_start = new_start + start
						new_end = new_end + start + 1

						# === next round: first line of the next gene ===
						this_start = max(astart, bstart)
						this_end = min(aend, bend)
						prev_cov_array = np.zeros(aend - astart)
						prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

						prev_astart = astart
						prev_aend = aend
						prev_geneid = ':'.join(map(str, [ageneid, astart, aend, achrom, astrand]))
						prev_chrom = achrom
						prev_strand = astrand
						prev_gene = this_gene

					count_genes_with_reads += 1
					if 'novel' not in geneid:
						count_genes_with_reads_annotated += 1

					if args.verbose:
						print '\ngene:', geneid, start, end, chrom, strand, new_start, new_end, \
							cov_array.size, str(datetime.now().time())

					if cov_array.size >= args.min_length:
						count_min_length_keep += 1

						# === get per-base coverage ===
						cov_avg = float(cov_sum) / float(cov_array.size)

						# === get introns ===
						intron_list = []
						exon_list = []
						jxn_list = []
						jxn_list, exon_list, intron_list = get_introns_exons(intron2jxnCount, cov_array.size, chrom, new_start, new_end, strand)

						if len(jxn_list) > 0:
							cov_avg_exon_with_utr, max_exon_cov = get_exon_cov(exon_list, cov_array)

							if strand != 0:
								strand_inferred = strand
							elif strand == 0:
								# === infer strand ===
								strand_inferred = infer_strand(intron_list, chrom, args.genome, args.verbose)

						else:
							cov_avg_exon_with_utr = cov_avg
							max_exon_cov = cov_avg
							strand_inferred = strand if strand != 0 else 'NA'

						if args.verbose:
							print '- junctions:', jxn_list
							print '- exons:', exon_list
							print '- strand:', strand_inferred
							print cov_avg_exon_with_utr, max_exon_cov

						if args.min_expn != -1:
							if cov_avg_exon_with_utr < args.min_expn:
								count_min_expn_filter += 1
								if 'novel' not in geneid:
									count_min_expn_filter_annotated += 1
								if args.verbose:
									print '- did not meet min expression level:', args.min_expn, cov_avg_exon_with_utr
								continue
							else:
								count_min_expn_keep += 1

						# === calculate cumulative vertical distance per gene ===
						vert_sum_array, vert_array = crs(cov_array)
						if vert_sum_array[0] == 'NA':
							count_vert_sum_max0 += 1
							continue

						# === KS test ===
						ksp = ks_test(vert_sum_array, args.plot, os.path.join(plot_dir, geneid))

						peak_inds_ttest_opt = []
						if ksp < args.test_thresh:
							vert2_array = vert_array * vert_array / max(vert_array)
							vert2_vert_sum_array, vert2_vert_array = crs(vert2_array)
							if vert2_vert_sum_array[0] == 'NA':
								count_vert2_sum_max0 += 1
								continue

							# === get parameters ===
							if args.winsize != -1:
								denoise_winsize_list = [args.winsize]
							else:
								winsize_max = max(100, int(round(vert_sum_array.size / 100, -2)))
								denoise_winsize_list = []
								winsize_min = 100
								while winsize_min <= min(winsize_max, 500):
									denoise_winsize_list.append(winsize_min)
									winsize_min = winsize_min + 100

								denoise_winsize_list = sorted(denoise_winsize_list)
							if args.verbose:
								print '- de-noising window sizes:', denoise_winsize_list

							amp_thresh_list = [args.peak_thresh] if args.peak_thresh != -1.0 else [0.05, 0.1, 0.15]
							if args.verbose:
								print '- amplitude thresholds:', amp_thresh_list

							peak_min_dist_list = [args.peak_min_dist] if args.peak_min_dist != -1 else [10, 50]
							if args.verbose:
								print '- min distance between peaks:', peak_min_dist_list

							h = 0
							njxns_detected = -100
							other_detected = 100
							param2cp = {}
							param2cpopt = {}
							param2fcopt = {}
							denoise_winsize_opt = 0
							amp_thresh_opt = 0
							peak_min_dist_opt = 0
							peak_inds_ttest = []
							peak_inds_ttest_with_ends = []
							for denoise_winsize in denoise_winsize_list:
								# get distance to line & denoise
								if args.verbose:
									print '- de-noising', denoise_winsize, str(datetime.now().time())
								line_dist_array2, line_dist_array2_denoise = get_linedist(vert2_vert_sum_array, denoise_winsize)
								line_dist_array, line_dist_array_denoise = get_linedist(vert_sum_array, denoise_winsize)

								a2totcp = {}
								for a, amp_thresh in enumerate(amp_thresh_list):
									if a > 0 and a2totcp[a - 1] == 0:
										if args.verbose:
											print '  - skipping higher amplitude thresholds because previous gave 0 change points'
										a2totcp[a] = 0
									else:
										for peak_min_dist in peak_min_dist_list:
											# --------------------------------------------------
											# peak calling
											# --------------------------------------------------
											if args.verbose:
												print '- parameters: window size', denoise_winsize, 'amplitude', amp_thresh, 'min distance', peak_min_dist, str(datetime.now().time())
											peak_inds_ip_combined_filtered = []
											peak_inds_ip_combined = []
											peak_inds2_ip_combined = []
											if np.unique(line_dist_array2_denoise).size != 1:
												# === crs^2 ===
												# call peaks
												peak_inds2 = peakutils.peak.indexes(line_dist_array2_denoise, thres=amp_thresh, min_dist=peak_min_dist)
												peak_inds_min2 = peakutils.peak.indexes(-1 * line_dist_array2_denoise, thres=amp_thresh, min_dist=peak_min_dist)

												# increase precision in case de-noising smoothed too much: call peaks on the real data just in the window around the peaks called on smooth data
												peak_inds2_ip = increase_precision(peak_inds2, denoise_winsize, line_dist_array2, 'max', amp_thresh)
												peak_inds2_min_ip = increase_precision(peak_inds_min2, denoise_winsize, line_dist_array2, 'min', amp_thresh)

												# combine
												peak_inds2_ip_combined = np.sort(np.append(peak_inds2_ip, peak_inds2_min_ip))
												if args.verbose:
													print '- peak inds max:', len(peak_inds2), str(datetime.now().time())
													print '- peak inds min:', len(peak_inds_min2), str(datetime.now().time())
													print '- increased precision max:', len(peak_inds2_ip), str(datetime.now().time())
													print '- increased precision min:', len(peak_inds2_min_ip), str(datetime.now().time())
													print '- peak inds crs^2:', len(peak_inds2_ip_combined), str(datetime.now().time())

											if np.unique(line_dist_array_denoise).size != 1:
												# === original crs ===
												# call peaks
												peak_inds = peakutils.peak.indexes(line_dist_array_denoise, thres=amp_thresh, min_dist=peak_min_dist)
												peak_inds_min = peakutils.peak.indexes(-1 * line_dist_array_denoise, thres=amp_thresh, min_dist=peak_min_dist)

												# increase precision in case de-noising smoothed too much: call peaks on the real data just in the window around the peaks called on smooth data
												peak_inds_ip = increase_precision(peak_inds, denoise_winsize, line_dist_array, 'max', amp_thresh)
												peak_inds_min_ip = increase_precision(peak_inds_min, denoise_winsize, line_dist_array, 'min', amp_thresh)

												# combine
												peak_inds_ip_combined = np.sort(np.append(peak_inds_ip, peak_inds_min_ip))
												if args.verbose:
													print '- peak inds max:', len(peak_inds), str(datetime.now().time())
													print '- peak inds min:', len(peak_inds_min), str(datetime.now().time())
													print '- increased precision max:', len(peak_inds_ip), str(datetime.now().time())
													print '- increased precision min:', len(peak_inds_min_ip), str(datetime.now().time())
													print '- peak inds crs:', len(peak_inds_ip_combined), str(datetime.now().time())

												if np.unique(line_dist_array2_denoise).size != 1:
													if peak_inds_ip_combined.size != 0 and peak_inds2_ip_combined.size != 0:
														peak_inds_ip_combined_filtered = np.append(peak_inds2_ip_combined, peak_inds_ip_combined[peak_inds_ip_combined < min(peak_inds2_ip_combined)])
														peak_inds_ip_combined_filtered = np.unique(np.append(peak_inds_ip_combined_filtered, peak_inds_ip_combined[peak_inds_ip_combined > max(peak_inds2_ip_combined)]))
													else:
														peak_inds_ip_combined_filtered = np.unique(np.append(peak_inds2_ip_combined, peak_inds_ip_combined))
												else:
													peak_inds_ip_combined_filtered = np.copy(peak_inds_ip_combined)

												# === UTRs: local crs ===
												if len(jxn_list) > 0:
													temp_cov_array = cov_array[:jxn_list[0]]
													if len(temp_cov_array) > 0:
														temp_vert_sum_array, temp_vert_array = crs(temp_cov_array)
														if temp_vert_sum_array[0] != 'NA' and len(temp_vert_sum_array) > 0 and max(temp_vert_sum_array) > 0:
															temp_ksp = ks_test(temp_vert_sum_array, args.plot, os.path.join(plot_dir, geneid) + '_temp1')
															if temp_ksp < args.test_thresh:
																temp_line_dist_array, temp_line_dist_array_denoise = get_linedist(temp_vert_sum_array, denoise_winsize)
																if np.unique(temp_line_dist_array_denoise).size != 1:
																	temp_peak_inds = peakutils.peak.indexes(temp_line_dist_array_denoise, thres=amp_thresh, min_dist=peak_min_dist)
																	temp_peak_inds_min = peakutils.peak.indexes(-1 * temp_line_dist_array_denoise, thres=amp_thresh, min_dist=peak_min_dist)
																	temp_peak_inds_ip = increase_precision(temp_peak_inds, denoise_winsize, temp_line_dist_array, 'max', amp_thresh)
																	temp_peak_inds_min_ip = increase_precision(temp_peak_inds_min, denoise_winsize, temp_line_dist_array, 'min', amp_thresh)
																	temp_peak_inds_ip_combined = np.append(temp_peak_inds_ip, temp_peak_inds_min_ip)
																	peak_inds_ip_combined_filtered = np.append(peak_inds_ip_combined_filtered, temp_peak_inds_ip_combined)

													temp_cov_array = cov_array[jxn_list[-1] - 1:]
													if len(temp_cov_array) > 0:
														temp_vert_sum_array, temp_vert_array = crs(temp_cov_array)
														if temp_vert_sum_array[0] != 'NA' and len(temp_vert_sum_array) > 0 and max(temp_vert_sum_array) > 0:
															temp_ksp = ks_test(temp_vert_sum_array, args.plot, os.path.join(plot_dir, geneid) + '_temp2')
															if temp_ksp < args.test_thresh:
																temp_line_dist_array, temp_line_dist_array_denoise = get_linedist(temp_vert_sum_array, denoise_winsize)
																if np.unique(temp_line_dist_array_denoise).size != 1:
																	temp_peak_inds = peakutils.peak.indexes(temp_line_dist_array_denoise, thres=amp_thresh, min_dist=peak_min_dist)
																	temp_peak_inds_min = peakutils.peak.indexes(-1 * temp_line_dist_array_denoise, thres=amp_thresh, min_dist=peak_min_dist)
																	temp_peak_inds_ip = increase_precision(temp_peak_inds, denoise_winsize, temp_line_dist_array, 'max', amp_thresh)
																	temp_peak_inds_min_ip = increase_precision(temp_peak_inds_min, denoise_winsize, temp_line_dist_array, 'min', amp_thresh)
																	temp_peak_inds_ip_combined = np.append(temp_peak_inds_ip, temp_peak_inds_min_ip)
																	peak_inds_ip_combined_filtered = np.append(peak_inds_ip_combined_filtered, temp_peak_inds_ip_combined)
											else:
												peak_inds_ip_combined_filtered = np.copy(peak_inds2_ip_combined)

											peak_inds_ip_combined_filtered = np.sort(np.unique(peak_inds_ip_combined_filtered))

											if args.verbose:
												print '- peak inds combined:', len(peak_inds_ip_combined_filtered), str(datetime.now().time())
												print '  ->', peak_inds_ip_combined_filtered

											# --------------------------------------------------
											# peak filtering
											# --------------------------------------------------
											if len(peak_inds_ip_combined_filtered) != 0:
												# === t-test: +/- window ===
												peak_inds_ttest, ind2tp = get_ttest(peak_inds_ip_combined_filtered, cov_array, denoise_winsize * 2, args.test_thresh)
												param2cp[(denoise_winsize, amp_thresh, peak_min_dist)] = peak_inds_ttest
												if args.verbose:
													print '- t-test:', len(peak_inds_ttest), str(datetime.now().time()), peak_inds_ttest

												# === filter #1: enforce peak_min_dist ===
												to_delete = []
												peak_inds_ttest_with_ends = [0] + peak_inds_ttest + [cov_array.size - 1]  # include ends
												for i, pos in enumerate(peak_inds_ttest_with_ends):
													if i > 0 and pos - peak_inds_ttest_with_ends[i - 1] <= peak_min_dist:
														ind = i - 1  # because we added the ends
														if i == 1:  # keep the end & delete the change point
															to_delete.append(ind)
														elif i == len(peak_inds_ttest_with_ends) - 1:  # keep the end & delete the change point
															to_delete.append(ind)
														else:  # keep the change point with lower t-test p-value
															to_delete.append(ind) if ind2tp[peak_inds_ttest[ind]] > ind2tp[peak_inds_ttest[ind - 1]] else to_delete.append(ind - 1)

												peak_inds_ttest = np.delete(peak_inds_ttest, to_delete)
												if args.verbose:
													print '- filtered by -d:', len(peak_inds_ttest), str(datetime.now().time()), peak_inds_ttest

												# === filter #2: fold change of first bps ===
												to_delete = []
												for i, peak in enumerate(peak_inds_ttest):
													if i > 0:
														prev_range_end = min(prev_peak + denoise_winsize, peak)
														this_range_end = min(peak + denoise_winsize, peak_inds_ttest[i + 1]) if i != len(peak_inds_ttest) - 1 else peak + denoise_winsize
														this_mean = np.mean(cov_array[peak:this_range_end])
														prev_mean = np.mean(cov_array[prev_peak:prev_range_end])
														log2fc = math.log((this_mean + 1) / (prev_mean + 1), 2)

														if abs(log2fc) < math.log(args.fcthresh, 2):
															if ind2tp[peak] > ind2tp[prev_peak]:
																to_delete.append(peak)
															else:
																to_delete.append(prev_peak)
																prev_peak = peak
														else:
															prev_peak = peak
													else:
														prev_peak = peak
												peak_inds_ttest = np.asarray([x for x in peak_inds_ttest if x not in to_delete])
												if args.verbose:
													print '- filtered by fold change of first bps:', len(peak_inds_ttest), str(datetime.now().time()), peak_inds_ttest

												# === filter #3: require increasing at 5' end & decreasing at 3' end ===
												if len(jxn_list) > 0:
													utr5 = peak_inds_ttest[peak_inds_ttest <= jxn_list[0] + args.juncdist]
													utr3 = peak_inds_ttest[peak_inds_ttest >= jxn_list[-1] - args.juncdist]
													peak_inds_ttest = get_end_cov(utr5, utr3, cov_array, peak_inds_ttest)
													if args.verbose:
														print '- filtered by 5\'/3\' ends:', len(peak_inds_ttest), str(datetime.now().time()), peak_inds_ttest

												# === filter #4: filter distal ends with low coverage ===
												to_delete = []
												peak_inds_ttest_with_ends = [0] + np.asarray(peak_inds_ttest).tolist() + [cov_array.size - 1]  # include ends
												if args.min_expn_distal != 0:
													if len(jxn_list) > 0:
														peak_inds_ttest_tandem_left = [x for x in peak_inds_ttest_with_ends if x <= int(jxn_list[0]) - args.juncdist]
														peak_inds_ttest_tandem_right = [x for x in peak_inds_ttest_with_ends if x >= int(jxn_list[-1]) + args.juncdist]

														# consider change points before the first junction
														if len(peak_inds_ttest_tandem_left) > 1:
															for i, ind in enumerate(peak_inds_ttest_tandem_left):
																next_cut = peak_inds_ttest_with_ends[i + 1]
																this_mean = round(np.mean(cov_array[ind:next_cut]), 2)
																if this_mean <= args.min_expn_distal:
																	to_delete.append(ind)

														# consider change points after the last junction
														if len(peak_inds_ttest_tandem_right) > 1:
															for i, ind in enumerate(peak_inds_ttest_tandem_right):
																prev_cut = peak_inds_ttest_with_ends[peak_inds_ttest_with_ends.index(ind) - 1]
																prev_mean = round(np.mean(cov_array[prev_cut:ind]), 2)
																if prev_mean <= args.min_expn_distal:
																	to_delete.append(ind)

														for ind in to_delete:
															peak_inds_ttest_with_ends.remove(ind)

													if args.verbose:
														print '- filter distal ends with low coverage:', len(peak_inds_ttest_with_ends), peak_inds_ttest_with_ends

												# === filter #5: fold change of whole segments ===
												to_delete = []
												fc_list = []
												if len(peak_inds_ttest_with_ends) > 0:
													for i, ind in enumerate(peak_inds_ttest_with_ends):
														next_cut = peak_inds_ttest_with_ends[i + 1] if i != len(peak_inds_ttest_with_ends) - 1 else cov_array.size
														if i == 0:
															this_mean = np.mean(cov_array[ind:next_cut])
															prev_mean = 0
															this_sd = np.std(cov_array[ind:next_cut])
															prev_sd = 0
														elif i == len(peak_inds_ttest_with_ends) - 1:
															this_mean = 0
															prev_mean = np.mean(cov_array[prev_cut:ind])
															this_sd = 0
															prev_sd = np.std(cov_array[prev_cut:ind])
														else:
															this_mean = np.mean(cov_array[ind:next_cut])
															prev_mean = np.mean(cov_array[prev_cut:ind])
															this_sd = np.std(cov_array[ind:next_cut])
															prev_sd = np.std(cov_array[prev_cut:ind])

														log2fc = math.log((this_mean + 1) / (prev_mean + 1), 2)

														if abs(log2fc) < math.log(args.fcthresh, 2) and i != 0 and i != len(peak_inds_ttest_with_ends) - 1:
															to_delete.append(i)
														else:
															prev_cut = ind

													peak_inds_ttest_with_ends = np.delete(peak_inds_ttest_with_ends, to_delete)

													# calculate fold change of whole segments
													for i, ind in enumerate(peak_inds_ttest_with_ends):
														next_cut = peak_inds_ttest_with_ends[i + 1] if i != len(peak_inds_ttest_with_ends) - 1 else cov_array.size
														if i == 0:
															this_mean = np.mean(cov_array[ind:next_cut])
															prev_mean = 0
															this_sd = np.std(cov_array[ind:next_cut])
															prev_sd = 0
														elif i == len(peak_inds_ttest_with_ends) - 1:
															this_mean = 0
															prev_mean = np.mean(cov_array[prev_cut:ind])
															this_sd = 0
															prev_sd = np.std(cov_array[prev_cut:ind])
														else:
															this_mean = np.mean(cov_array[ind:next_cut])
															prev_mean = np.mean(cov_array[prev_cut:ind])
															this_sd = np.std(cov_array[ind:next_cut])
															prev_sd = np.std(cov_array[prev_cut:ind])

														log2fc = math.log((this_mean + 1) / (prev_mean + 1), 2)
														prev_cut = ind
														fc_list.append(':'.join(map(str, [this_mean, log2fc])))
												else:
													count_filter123 += 1

												if args.verbose:
													print '- fold change whole segment:', len(peak_inds_ttest_with_ends), peak_inds_ttest_with_ends

												param2cpopt[(denoise_winsize, amp_thresh, peak_min_dist)] = peak_inds_ttest_with_ends
												param2fcopt[(denoise_winsize, amp_thresh, peak_min_dist)] = fc_list
											elif args.verbose:
												count_no_cps += 1
												print '  - no change points called'

											# get total change points for this amp threshold
											if a not in a2totcp:
												a2totcp[a] = len(peak_inds_ttest_with_ends)
											elif len(peak_inds_ttest_with_ends) > a2totcp[a]:
												a2totcp[a] = len(peak_inds_ttest_with_ends)

											# optimization criteria
											if len(peak_inds_ttest_with_ends) > 0:
												if len(jxn_list) > 0:
													# optimal parameters: detect junctions within 50bp
													jxn_closest = [min(peak_inds_ttest_with_ends, key=lambda x:abs(x - j)) for j in jxn_list]
													jxn_peak_dif = [abs(j - c) for j, c in zip(jxn_list, jxn_closest)]
													tot_near_jxns = len([x for x in jxn_peak_dif if x <= args.juncdist])
													tot_other = len(peak_inds_ttest_with_ends) - tot_near_jxns

													if tot_near_jxns - tot_other > njxns_detected:
														njxns_detected = tot_near_jxns - tot_other
														denoise_winsize_opt = denoise_winsize
														amp_thresh_opt = amp_thresh
														peak_min_dist_opt = peak_min_dist
													elif args.winsize != -1:
														njxns_detected = tot_near_jxns - tot_other
														denoise_winsize_opt = denoise_winsize
														amp_thresh_opt = amp_thresh
														peak_min_dist_opt = peak_min_dist
												else:
													# optimal parameters: fewest total change points
													if len(peak_inds_ttest_with_ends) < other_detected:
														other_detected = len(peak_inds_ttest_with_ends)
														denoise_winsize_opt = denoise_winsize
														amp_thresh_opt = amp_thresh
														peak_min_dist_opt = peak_min_dist
											else:
												count_no_cps_ttest0 += 1
												if args.verbose:
													print '  - no change points found: total passing t-test = 0'

							# --------------------------------------------------
							# optimal change points
							# --------------------------------------------------
							if denoise_winsize_opt == 0 or amp_thresh_opt == 0 or peak_min_dist_opt == 0:
								count_no_cps_afterfiltering += 1
								if args.verbose:
									print '- no change points called', str(datetime.now().time())
							else:
								count_cps_called += 1
								if args.verbose:
									print '- optimal parameters: jxns detected', njxns_detected, '/', len(jxn_list), ': window size', denoise_winsize_opt, 'amplitude', amp_thresh_opt, 'min distance', peak_min_dist_opt, str(datetime.now().time())

								peak_inds_ttest_opt = param2cpopt[(denoise_winsize_opt, amp_thresh_opt, peak_min_dist_opt)].tolist()
								peak_inds_ttest_preopt = param2cp[(denoise_winsize_opt, amp_thresh_opt, peak_min_dist_opt)]
								fc_list_opt = param2fcopt[(denoise_winsize_opt, amp_thresh_opt, peak_min_dist_opt)]

								# === filter low coverage ends ===
								if args.max_end_ru != 0:
									cov_mean_list = [float(x.split(':')[0]) for x in fc_list_opt]
									ru_all = [x / max(cov_mean_list) for x in cov_mean_list]

									while ru_all[0] < args.max_end_ru:
										if args.verbose:
											print '-> removing left end:', ru_all[0], peak_inds_ttest_opt[0] + new_start
										del ru_all[0]
										del peak_inds_ttest_opt[0]
										del fc_list_opt[0]

									if len(ru_all) > 1:
										while ru_all[-2] < args.max_end_ru:
											if args.verbose:
												print '-> removing right end:', ru_all[-1], peak_inds_ttest_opt[-1] + new_start
											del ru_all[-1]
											del peak_inds_ttest_opt[-1]
											del fc_list_opt[-1]

								# label change points
								for i, peak_ind in enumerate(peak_inds_ttest_opt):
									if i == 0:
										if strand_inferred == '+':
											label = 'DistalTSS'
										elif strand_inferred == '-':
											label = 'DistalPolyA'
										else:
											label = 'DistalLeft'
									elif i == len(peak_inds_ttest_opt) - 1:
										if strand_inferred == '+':
											label = 'DistalPolyA'
										elif strand_inferred == '-':
											label = 'DistalTSS'
										else:
											label = 'DistalRight'
									elif len(jxn_list) > 0:
										jxn_distances = [abs(int(peak_ind) - int(x)) for x in jxn_list]
										exon_start_end_list = [map(int, exon.split(':')) for exon in exon_list]
										exon_flag_list = [1 if estart <= int(peak_ind) < eend else 0 for estart, eend in exon_start_end_list]

										if any(x <= args.juncdist for x in jxn_distances):
											label = 'Junction'
										elif int(peak_ind) <= int(jxn_list[0]) - args.juncdist:
											if strand_inferred == '-':
												label = 'TandemAPA'
											elif strand_inferred == '+':
												label = 'TandemATSS'
											else:
												label = 'TandemLeft'
										elif int(peak_ind) >= int(jxn_list[-1]) + args.juncdist:
											if strand_inferred == '+':
												label = 'TandemAPA'
											elif strand_inferred == '-':
												label = 'TandemATSS'
											else:
												label = 'TandemRight'
										else:
											if sum(exon_flag_list) == 0:
												label = 'Intron'
											else:
												label = 'Exon'
									else:
										label = 'Exon'

									if int(peak_ind) > cov_array.size - 1:
										print '\ngene:', geneid, start, end, chrom, strand, new_start, new_end, cov_array.size, str(datetime.now().time())
										print cov_array.size, peak_ind
										sys.exit(1)

									this_gs = geneid.split(':')[0]
									this_pos = str(int(peak_ind) + new_start)
									this_fc = fc_list_opt[i]
									cov_exons, fc = this_fc.split(':')

									# === write output ===
									if strand == '1':
										strand = '+'
									if strand == '-1':
										strand = '-'
									strand_inferred = strand if strand == '+' or strand == '-' else strand_inferred

									name = ':'.join(map(str, [label, this_gs, start, end, strand_inferred,
															  denoise_winsize_opt,
															  cov_exons,
															  cov_avg_exon_with_utr]))

									# the change point is 1-based before the change -> adjust by 1bp if change point goes from high to low so that reported change points are outside of the higher-coverage region
									if float(fc) < 0:
										this_pos = int(this_pos) - 1

									if strand != '+' and strand != '-':
										if strand_inferred == '+' or strand_inferred == '-':
											o.write('\t'.join([chrom, str(int(this_pos)), str(int(this_pos) + 1), name, fc, strand_inferred]) + '\n')
										else:
											o.write('\t'.join([chrom, str(int(this_pos)), str(int(this_pos) + 1), name, fc]) + '\n')
									else:
										o.write('\t'.join([chrom, str(int(this_pos)), str(int(this_pos) + 1), name, fc, strand]) + '\n')
						else:
							count_high_ksp += 1
							if args.verbose:
								print '- no change points called: KS p =', ksp

						# --------------------------------------------------
						# plots
						# --------------------------------------------------
						if args.plot:
							if args.verbose:
								print '- plotting'
							# plot coverage without change points
							x = np.linspace(0, len(cov_list) - 1, len(cov_list)) / 1000
							y = np.asarray(cov_list)
							out_plot = os.path.join(plot_dir, geneid + '_cov.pdf')
							pdf = PdfPages(out_plot)
							fig = pyplot.figure(figsize=(3, 3))
							pyplot.plot(x, y, color='k')
							pyplot.xlabel('Position (kb)')
							pyplot.ylabel('Coverage')
							pyplot.gcf().subplots_adjust(bottom=0.3, left=0.3)
							pdf.savefig()
							pdf.close()
							pyplot.close(fig)

							if len(jxn_list) > 0:
								# plot coverage without change points, excluding introns
								x = np.linspace(0, len(temp_cov_list) - 1, len(temp_cov_list)) / 1000
								y = np.asarray(temp_cov_list)
								out_plot = os.path.join(plot_dir, geneid + '_cov_noIntrons.pdf')
								pdf = PdfPages(out_plot)
								fig = pyplot.figure(figsize=(3, 3))
								pyplot.plot(x, y, color='k')
								pyplot.xlabel('Position (kb)')
								pyplot.ylabel('Coverage')
								pyplot.gcf().subplots_adjust(bottom=0.3, left=0.3)
								pdf.savefig()
								pdf.close()
								pyplot.close(fig)

							# plot coverage (black) + + change points (red)
							x = np.linspace(0, len(cov_list) - 1, len(cov_list)) / 1000
							y = np.asarray(cov_list)
							out_plot = os.path.join(plot_dir, geneid + '_cov_segs.pdf')
							pdf = PdfPages(out_plot)
							fig = pyplot.figure(figsize=(3, 3))
							pyplot.plot(x, y, color='k')
							if len(peak_inds_ttest_opt) > 0:
								for seg in peak_inds_ttest_opt:
									seg = float(seg) / float(1000)
									pyplot.axvline(x=seg, color='r')
							pyplot.title('Segmentation: n = ' + str(len(peak_inds_ttest_opt)))
							pyplot.xlabel('Position (kb)')
							pyplot.ylabel('Coverage')
							pyplot.gcf().subplots_adjust(bottom=0.3, left=0.3)
							pdf.savefig()
							pdf.close()
							pyplot.close(fig)

							if ksp < args.test_thresh:
								x = np.linspace(0, len(line_dist_array_denoise) - 1, len(line_dist_array_denoise)) / 1000
								y1 = line_dist_array2_denoise
								y2 = line_dist_array2
								y3 = line_dist_array_denoise
								y4 = line_dist_array
								out_plot = os.path.join(plot_dir, geneid + '_crsToLine.pdf')
								pdf = PdfPages(out_plot)
								fig = pyplot.figure(figsize=(3, 3))
								pyplot.plot(x, y1, color='b')
								pyplot.plot(x, y2, color='g')
								pyplot.plot(x, y3, color='c')
								pyplot.plot(x, y4, color='k')
								pyplot.title('crs to line: n = ' + str(len(peak_inds_ttest_preopt)))
								pyplot.xlabel('Position (kb)')
								pyplot.ylabel('Distance to line y=ax')
								pyplot.gcf().subplots_adjust(bottom=0.3, left=0.3)
								pdf.savefig()
								pdf.close()
								pyplot.close(fig)

								x = np.linspace(0, len(line_dist_array_denoise) - 1, len(line_dist_array_denoise)) / 1000
								y1 = line_dist_array2_denoise
								y2 = line_dist_array2
								y3 = line_dist_array_denoise
								y4 = line_dist_array
								out_plot = os.path.join(plot_dir, geneid + '_crsToLine_segs.pdf')
								pdf = PdfPages(out_plot)
								fig = pyplot.figure(figsize=(3, 3))
								pyplot.plot(x, y1, color='b')
								pyplot.plot(x, y2, color='g')
								pyplot.plot(x, y3, color='c')
								pyplot.plot(x, y4, color='k')
								if len(peak_inds_ttest_preopt) > 0:
									for seg in peak_inds_ttest_preopt:
										seg = float(seg) / float(1000)
										pyplot.axvline(x=seg, color='r')
								pyplot.title('crs to line: n = ' + str(len(peak_inds_ttest_preopt)))
								pyplot.xlabel('Position (kb)')
								pyplot.ylabel('Distance to line y=ax')
								pyplot.gcf().subplots_adjust(bottom=0.3, left=0.3)
								pdf.savefig()
								pdf.close()
								pyplot.close(fig)

					else:
						count_min_length_filter += 1
						if 'novel' not in geneid:
							count_min_length_filter_annotated += 1
						if args.verbose:
							print 'gene length <', args.min_length, '-> skipping'

	o.close()


	sort_bedfile(args.output + '.temp', args.output)
	os.remove(args.output + '.temp')

	if args.verbose:
		print count_genes_with_reads, 'total TUs with reads,', count_genes_with_reads_annotated, 'annotated'
		print count_min_length_filter, 'total TUs filtered by min TU length >=', args.min_length, '(', count_min_length_filter_annotated, 'annotated),', count_min_length_keep, 'kept'
		print count_min_expn_filter, 'total TUs filtered by min expression >=', args.min_expn, '(', count_min_expn_filter_annotated, 'annotated), kept', count_min_expn_keep
		print count_high_ksp, 'total TUs with ksp >=', args.test_thresh
		print count_vert_sum_max0, 'total TUs filtered by having max CVS = 0'
		print count_vert2_sum_max0, 'total TUs filtered by having max CVS2 = 0'
		print count_no_cps_ttest0, 'total TUs filtered by having no cps passing ttest'
		print count_filter123, 'total filtered by filters 1-3'
		print count_no_cps_afterfiltering, 'total TUs with no cps after all filters'
		print count_cps_called, 'total TUs with change points'

	# === delete temporary files ===
	os.remove(bgfile)
	print '\nfinished:', str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
