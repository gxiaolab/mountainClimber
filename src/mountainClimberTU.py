#!/user/bin/python -tt
"""
Call transcription units (TUs) de novo from RNA-Seq in each individual sample.
There are three major steps:
1. keep bins of window size w with >= n average read density (count per bp / length)
2. merge consecutive bins if >= p% of positions covered
3. trim/extend ends
"""


import os
import sys
import argparse
import gzip
import re
import math
import pybedtools as pb
from datetime import datetime
from collections import defaultdict
from functions import run_command, sort_bedfile


def scientific2int(count):
	"""Convert scientific notation to integer"""
	(num, exp) = count.split('e+')
	count = int(float(num) * float(math.pow(10, int(exp))))
	return count


def window_count(bedgraph, window_size):
	"""Count the breadth, depth, and length of each window"""
	window2count_length = {}
	chr2starts = defaultdict(list)
	chr2ends = defaultdict(list)

	f = gzip.open(bedgraph, 'rb') if re.search(r'\.gz$', bedgraph) else open(bedgraph, 'r')

	# get coverage & length for each window
	last_chrom = -1
	last_start = -1
	last_end = -1
	for line in f:
		if not line.startswith('track') and not line.startswith('browser'):
			x = line.rstrip().split('\t')
			if len(x) != 4:
				sys.stderr.write('EXIT: do not recognize format! Please input bedgraph (4 fields)\n')
				sys.stderr.write(x)
				sys.exit(1)
			(chrom, start, end, count) = x
			start = int(start)
			end = int(end)
			count = scientific2int(count) if '+' in count else int(count)

			if count != 0:  # skip regions with zero coverage
				# === get continuous regions of coverage ===
				if last_start == -1:  # first start
					chr2starts[chrom].append(start)
				elif prev_chrom == chrom:
					if last_end != start:  # new block in same chrom
						chr2starts[chrom].append(start)
						chr2ends[chrom].append(last_end)
				else:  # new chrom
					chr2starts[chrom].append(start)
					chr2ends[prev_chrom].append(last_end)
				prev_chrom = chrom
				last_start = start
				last_end = end

				# === get 1kb windows ===
				window_lower = int(round(start / window_size) * window_size)
				window_upper = end + (window_size - end) % window_size
				total_windows = (window_upper - window_lower) / window_size

				if (total_windows == 1):  # region spans one window
					length = end - start
					total_count = count * length

					# update window coverage & length
					if (chrom, window_lower, window_upper) not in window2count_length:
						window2count_length[(chrom, window_lower, window_upper)] = (total_count, length, start, end)
					else:
						(prev_count, prev_length, prev_start, prev_end) = window2count_length[(chrom, window_lower, window_upper)]
						window2count_length[(chrom, window_lower, window_upper)] = (prev_count + total_count, prev_length + length, min(prev_start, start), max(prev_end, end))
				else:  # region spans multiple windows
					for window_start in xrange(window_lower, window_upper, window_size):
						window_end = window_start + window_size

						# get coverage
						if window_start <= start:
							# first window
							this_start = start
							this_end = window_end
						elif window_end == window_upper:
							# last window
							this_start = window_start
							this_end = end
						else:
							# middle windows
							this_start = window_start
							this_end = window_end

						# update window coverage & length
						length = this_end - this_start
						total_count = count * length

						if (chrom, window_start, window_end) not in window2count_length:
							window2count_length[(chrom, window_start, window_end)] = (total_count, length, this_start, this_end)
						else:
							(prev_count, prev_length, prev_start, prev_end) = window2count_length[(chrom, window_start, window_end)]
							window2count_length[(chrom, window_start, window_end)] = (prev_count + total_count, prev_length + length, min(prev_start, this_start), max(prev_end, this_end))
	f.close()
	chr2ends[chrom].append(end)  # last end
	return window2count_length, chr2starts, chr2ends


def merge_windows(window2count_length, min_percent, min_reads, bedgraph, outfile_bgmerge, outfile_temp_windows):
	"""
	Merge consecutive windows that meet the following criteria:
			>= min_reads average reads/bp in the window
			>= min_percent bases in the window covered by reads
	Report window start as the first nucleotide with reads (which may be before or after window start)
	Report window end as the first nucleotide with reads (which may be before or after window end)
	"""
	transcripts_trimmed = {}

	# === filter windows & merge ===
	count_windows_filtered = 0
	windows = list(sorted(window2count_length.keys()))
	t = open(outfile_temp_windows, 'w')
	for window in windows:
		(count, length, rstart, rend) = window2count_length[window]
		(chrom, wstart, wend) = window
		p = float(length) / (float(wend) - float(wstart))
		n = float(count) / float(length)

		if n >= min_reads and p >= min_percent:
			count_windows_filtered += 1
			t.write('\t'.join(map(str, [chrom, wstart, wend, ':'.join(map(str, [rstart, rend, count, length]))])) + '\n')
	t.close()

	if count_windows_filtered == 0:
		sys.stderr.write(' '.join(map(str, ['  -> No windows met criteria: >=', min_reads, 'reads/bp and >=', min_percent, "% bp covered"])) + '\n')
		return 0
	else:
		outfile_merged = outfile_temp_windows + '.merge'
		pb.BedTool(outfile_temp_windows).sort().merge(c=4, o='distinct').saveas(outfile_merged)

		# merge bedgraph
		pb.BedTool(bedgraph).merge().saveas(outfile_bgmerge)

		# merge TUs with bedgraph
		outfile_extend = outfile_merged + '.extend'
		pb.BedTool(outfile_merged).sort().intersect(wao=True, b=outfile_bgmerge).saveas(outfile_extend)

		tx2start = {}
		tx2end = {}
		with open(outfile_extend, 'r') as f:
			for line in f:
				(tchrom, tstart, tend, tname, bchrom, bstart, bend, bcov) = line.rstrip().split('\t')
				bstart, bend = map(int, [bstart, bend])

				if (tchrom, tstart, tend) not in tx2start:
					tx2start[(tchrom, tstart, tend)] = bstart
				elif bstart < tx2start[(tchrom, tstart, tend)]:
					tx2start[(tchrom, tstart, tend)] = bstart

				if (tchrom, tstart, tend) not in tx2end:
					tx2end[(tchrom, tstart, tend)] = bend
				elif bend > tx2end[(tchrom, tstart, tend)]:
					tx2end[(tchrom, tstart, tend)] = bend

		for (tchrom, tstart, tend) in tx2start:
			new_start = tx2start[(tchrom, tstart, tend)]
			new_end = tx2end[(tchrom, tstart, tend)]
			transcripts_trimmed[(tchrom, new_start, new_end)] = 1

		# --- cleanup ---
		os.remove(outfile_extend)
		os.remove(outfile_bgmerge)
		os.remove(outfile_merged)

		return transcripts_trimmed


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									 description='Call transcription units (TUs) de novo from RNA-Seq.')
	group = parser.add_argument_group('Input')
	group.add_argument('-b', '--bedgraph', dest='bedgraph', type=str, metavar='',
					   help='Bedgraph file. Can be .gz compressed.')
	group.add_argument('-j', '--junc', dest='junc', type=str, metavar='',
					   help='Junction .bedgraph or .bed file. If file suffix is .bed, will convert to bedgraph.')
	group.add_argument('-g', '--genome', dest='genome', type=str, metavar='',
					   help='Input chromosome sizes.')

	group = parser.add_argument_group('Parameters')
	group.add_argument('-c', '--minjxncount', dest='minjxncount', type=int, default=2, metavar='',
					   help='Minimum junction read count.')
	group.add_argument('-s', '--strand', dest='strand', type=int, choices=[1, -1, 0], required=True,
					   help='Strand of bedgraph file.')
	group.add_argument('-w', '--window_size', dest='window_size', type=int, default=1000, metavar='',
					   help='Window size.')
	group.add_argument('-p', '--min_percent', dest='min_percent', type=float, default=1.0, metavar='',
					   help='Minimum percentage of the window covered.')
	group.add_argument('-n', '--min_reads', dest='min_reads', type=int, default=10, metavar='',
					   help='Minimum number of reads per window.')

	group = parser.add_argument_group('Output')
	group.add_argument('-o', '--output', dest='output', type=str, metavar='',
					   help='Output bed filename.')
	args = parser.parse_args()

	# --------------------------------------------------
	# main routine
	# --------------------------------------------------
	print '\njob starting:', str(datetime.now().time())

	if not args.output:
		sys.stderr.write('EXIT: please enter --output')
		sys.exit(1)

	if not args.bedgraph:
		sys.stderr.write('EXIT: please enter --bedgraph')
		sys.exit(1)

	# === calculate window read counts. incorporate introns if --junc provided. ===
	if args.junc:
		if re.search(r'\.bed$', args.junc):
			print '- converting junction .bed file to .bedgraph', str(datetime.now().time())

			# --- get overlapping junction regions ---
			temp_jxn_intersect = args.output + '.temp.junction.intersect.txt'
			if args.strand == 0:
				temp_jxn_file = args.junc
				pb.BedTool(args.junc).genome_coverage(bg=True, g=args.genome).intersect(
					wo=True, b=temp_jxn_file).saveas(temp_jxn_intersect)
			else:
				if args.strand == 1:
					gc_strand = '+'
				elif args.strand == -1:
					gc_strand = '-'

				# get this strand's junctions only & filter junctions by read counts
				temp_jxn_file = args.output + '.' + str(args.strand) + '.jxn.bed'
				o = open(temp_jxn_file, 'w')
				with open(args.junc, 'r') as f:
					for line in f:
						(chrom, start, end, name, cov, strand) = line.rstrip().split('\t')
						if strand == gc_strand and int(cov) >= args.minjxncount:
							o.write(line)
				o.close()

				# get junction genomecov & intersect to get the # reads per genomecov region
				pb.BedTool(args.junc).genome_coverage(bg=True, g=args.genome, strand=gc_strand).intersect(wo=True, b=temp_jxn_file).saveas(temp_jxn_intersect)

			# --- get bedgraph of junction regions ---
			jxn2cov = {}
			with open(temp_jxn_intersect, 'r') as f:
				for line in f:
					if args.strand == 0:
						(cchrom, cstart, cend, ccov, jchrom, jstart, jend, jname, jcov, overlap) = line.rstrip().split('\t')
					else:
						(cchrom, cstart, cend, ccov, jchrom, jstart, jend, jname, jcov, jstrand, overlap) = line.rstrip().split('\t')

					jxn = ':'.join([cchrom, cstart, cend])
					if jxn not in jxn2cov:
						jxn2cov[jxn] = int(jcov)
					else:
						jxn2cov[jxn] += int(jcov)

			outfile_temp_jxn_bg_unsorted = args.output + '.' + str(args.strand) + '.jxn.bedgraph.unsorted'
			o = open(outfile_temp_jxn_bg_unsorted, 'w')
			for jxn in jxn2cov:
				o.write('\t'.join(jxn.split(':')) + '\t' + str(jxn2cov[jxn]) + '\n')
			o.close()

			outfile_temp_jxn_bg = args.output + '.' + str(args.strand) + '.jxn.bedgraph'
			sort_bedfile(outfile_temp_jxn_bg_unsorted, outfile_temp_jxn_bg)

			# --- cleanup ---
			os.remove(temp_jxn_intersect)
			os.remove(outfile_temp_jxn_bg_unsorted)
			if args.strand != 0:
				os.remove(temp_jxn_file)

		elif re.search(r'\.bedgraph$', args.junc):
			outfile_temp_jxn_bg = args.junc
		else:
			sys.stderr.write('EXIT: File extension not recognized. Please input .bed or .bedgraph --junc file')
			sys.exit(1)

		# combine bedgraphs: using subprocess because pybedtools complains about unionbedg format not being bed
		print '- merging junctions with bedgraph', str(datetime.now().time())
		outfile_temp_unionbg = args.output + '.unionbedg'
		j = open(outfile_temp_unionbg, 'w')
		cmd = ['bedtools', 'unionbedg', '-i', args.bedgraph, outfile_temp_jxn_bg]
		stdout, stderr = run_command(cmd, stdoutfile=j)
		j.close()
		os.remove(outfile_temp_jxn_bg)

		# write output
		outfile_temp_merge_bg_jxn = args.output + '.bg_jxn.bedgraph'
		m = open(outfile_temp_merge_bg_jxn, 'w')
		with open(outfile_temp_unionbg, 'r') as f:
			for line in f:
				(chrom, start, end, count1, count2) = line.split('\t')
				if '+' in count1 or '+' in count2:
					count1 = scientific2int(count1) if '+' in count1 else int(count1)
					count2 = scientific2int(count2) if '+' in count2 else int(count2)
				m.write('\t'.join(map(str, [chrom, start, end, int(count1) + int(count2)])) + '\n')
		m.close()

		print '- calculating window read counts:', str(datetime.now().time())
		window2count_length, chr2starts, chr2ends = window_count(outfile_temp_merge_bg_jxn, args.window_size)
	else:
		sys.stderr.write('No --junc file was input. It is recommended to include this input.')

		print '- calculating window read counts:', str(datetime.now().time())
		window2count_length, chr2starts, chr2ends = window_count(args.bedgraph, args.window_size)
		outfile_temp_merge_bg_jxn = args.bedgraph

	print '  ->', len(window2count_length.keys()), 'total windows', str(datetime.now().time())

	# === merge consecutive windows ===
	print '- merging consecutive windows:', str(datetime.now().time())
	outfile_temp_windows = args.output + '.windows'
	outfile_temp_bg_merged = args.output + '.bg.merged'
	transcripts_trimmed = merge_windows(window2count_length, args.min_percent, args.min_reads, outfile_temp_merge_bg_jxn, outfile_temp_bg_merged, outfile_temp_windows)

	if transcripts_trimmed != 0:
		print '  ->', len(transcripts_trimmed.keys()), 'total transcripts after merging'

		# === get strand ===
		if args.strand == 1:
			strand = "+"
		elif args.strand == -1:
			strand = "-"
		elif args.strand == 0:
			strand = 0
		else:
			print 'DIED: did not recognize strand', args.strand
			sys.exit()

		# === print output -> unsorted bed format ===
		print "printing output (unsorted): " + str(datetime.now().time())
		transcripts_trimmed_list = list(sorted(transcripts_trimmed.keys()))
		output_temp = args.output + '.temp'
		o = open(output_temp, 'w')
		for ind, tx in enumerate(transcripts_trimmed_list):
			(chrom, start, end) = tx
			if not chrom.startswith('chr'):
				chrom = 'chr' + chrom
			if strand != 0:
				o.write('\t'.join([chrom, str(start), str(end), 'tx_' + str(ind), '0', strand]) + '\n')
			else:
				o.write('\t'.join([chrom, str(start), str(end), 'tx_' + str(ind)]) + '\n')
		o.close()

		# === output: sorted bed format & delete the temp unsorted file ===
		sort_bedfile(output_temp, args.output)

		# === clean up ===
		os.remove(output_temp)

	# === clean up ===
	os.remove(outfile_temp_windows)
	if args.junc:
		os.remove(outfile_temp_unionbg)
		os.remove(outfile_temp_merge_bg_jxn)

	print '\nfinished: ' + str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
