#!/user/bin/python -tt
"""
Identify split reads from a bam file and report the total counts per intron.
"""


import os
import sys
import argparse
from datetime import datetime
import pysam
from collections import OrderedDict


def make_intron_bed(bam_file, outfile, overhang, min_intron, max_intron, strandedness):
	"""
	Input: bam_file
	Procedure:
		Identify reads with gaps in alignment.
		Report those with:
			gap length between min_intron and max_intron
			>= overhang nucleotides matched
	Output: outfile in bed format: chromosome, intron start, intron end, intron name, read count, strand
	"""

	bamfile = pysam.Samfile(bam_file, "rb")
	out = open(outfile + '.tmp','w')
	overhang, min_intron, max_intron = map(int, (overhang, min_intron, max_intron))

	if strandedness == 'fr-firststrand':
		strnd1 = '-' ; strnd2 = '+'
	elif strandedness == 'fr-secondstrand':
		strnd1 = '+' ; strnd2 = '-'
	elif strandedness == 'fr-unstrand':
		strnd1 = 'nss' ; strnd2 = 'nss'
	elif strandedness == 'single':
		strnd1 = '+' ; strnd2 = '-'
	else:
		sys.stderr.write('EXIT: do not recongize strandedness ' + strandedness + '\n')
		sys.exit(1)

	junctions_idme = OrderedDict()
	for read in bamfile:
		if len(read.cigar) > 0 and 'N' in read.cigarstring:
			chrm = read.reference_name
			if strandedness == 'single':
				if read.is_reverse: strand = strnd1
				else: strand = strnd2
			else:
				if read.is_read1 and not read.is_reverse: strand = strnd1
				elif read.is_read2 and read.is_reverse: strand = strnd1
				else: strand = strnd2

			match_lengths, block_lengths, junct_lengths = process_CIGAR(read.cigar)

			for i in range(len(block_lengths)-1):
				junc_stt = read.pos + block_lengths[i]		# 0-based, first pos on intron
				junc_end = junc_stt + junct_lengths[i]	 	# 1-based, last pos on intron
				junc_len = junct_lengths[i]
				junction = ':'.join([chrm, str(junc_stt), str(junc_end), strand])
				if min_intron <= junc_len <= max_intron:
					if match_lengths[i] >= overhang and match_lengths[i+1] >= overhang:
						try: 	junctions_idme[junction][4] += 1
						except:	junctions_idme[junction] = [chrm, junc_stt, junc_end, junction, 1, strand]

	for line in junctions_idme.values():
		outline = map(str, line)
		if strnd1 == 'nss':
			out.write('\t'.join(outline[:-1]) + '\n')
		else:
			out.write('\t'.join(outline) + '\n')
	bamfile.close()
	out.close()

	os.system('cat %s.tmp | sort -k1,1 -k2,2n > %s' %(outfile, outfile))
	os.system('rm  %s.tmp' %(outfile))
	# To run separate chromosomes at the time or all at once
	chrm_number = list(set(zip(*junctions_idme.values())[0]))

	if len(chrm_number) == 1 : return chrm_number[0], strnd1
	elif len(chrm_number) > 1: return 'all', strnd1


def process_CIGAR(cigar):
	"""
	Process the sam file CIGAR string:
	0 => matches   (M)
	1 => insertion (I)
	2 => deletion  (D)
	3 => intron    (N)
	4 => soft clip (S)
	5 => hard clip (H)
	"""
	match_lengths = []
	block_lengths = []
	junct_lengths = []
	matches = 0 ; blocks = 0 ; no_of_blocks = 1
	for (ctype, length) in cigar:
		if ctype == 0:
			matches += length
			blocks  += length
		elif ctype == 3:
			junct_lengths.append(length)
			match_lengths.append(matches)
			block_lengths.append(blocks)
			blocks += length
			matches = 0
			no_of_blocks += 1
	block_lengths.append(blocks)
	match_lengths.append(matches)
	assert len(block_lengths) == no_of_blocks
	assert len(match_lengths) == no_of_blocks
	return match_lengths, block_lengths, junct_lengths


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	strand_choice_list = ['fr-firststrand', 'fr-secondstrand', 'fr-unstrand', 'single']

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Identify split reads from a bam file and report the total counts per intron.')
	group = parser.add_argument_group('Input')
	group.add_argument('-i', '--input_bam', dest='input_bam', type=str, metavar='', required=True,
		help='Bam file')
	group.add_argument('-s', '--strand', dest='strand', type=str, default = 'fr-firststrand',
		choices=strand_choice_list, metavar='', required=True,
		help='Strandedness. Options: ' + ', '.join(strand_choice_list))

	group = parser.add_argument_group('Parameters')
	group.add_argument('-e', '--overhang', dest='overhang', type=int, default=8, metavar='',
		help='Minimum number of base pairs in each exon')
	group.add_argument('-m', '--min_intron', dest='min_intron', type=int, default=30, metavar='',
		help='Minimum intron length')
	group.add_argument('-a', '--max_intron', dest='max_intron', type=int, default=500000, metavar='',
		help='Maximum intron length')

	group = parser.add_argument_group('Output')
	group.add_argument('-o', '--output', dest='output', type=str, required=True, metavar='',
		help='Output filename')
	args = parser.parse_args()

	# --------------------------------------------------
	# main routine
	# --------------------------------------------------
	print '\nstarting:', str(datetime.now().time())

	make_intron_bed(args.input_bam, args.output, overhang = args.overhang,
		min_intron = args.min_intron, max_intron = args.max_intron, strandedness = args.strand)

	print '\nfinished:', str(datetime.now().time())

# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
