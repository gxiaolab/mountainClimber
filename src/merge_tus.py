#!/user/bin/python -tt
"""
Merge TUs that overlap across samples and merge those within the same gene annotation.
"""

import os
import re
import sys
import argparse
import pybedtools as pb
from datetime import datetime
from collections import defaultdict


def check_empty_file(bedfile):
	"""Check if file is empty"""
	f = open(bedfile, 'r')
	if len(f.readlines()) == 0:
		sys.stderr.write('EXIT: no bedgraph output!')
		sys.exit(1)
	f.close()


def merge_annotation(refgtf, ss, temp_annotbed, temp_annotbed_sort, temp_annotbed_merged,
					 temp_annotbed_merged_sorted):
	"""Merge overlapping gene annotations & create duplicate gene IDs if genes are duplicated"""
	gs2minTxStart = {}
	gs2maxTxEnd = {}
	with open(refgtf, 'r') as g:
		for line in g:
			if line.startswith('#'):
				continue
			x = line.rstrip().split('\t')
			if x[2] == 'transcript':
				gs = [a for a in x[8].split(';') if 'gene_name' in a]
				tx = [a for a in x[8].split(';') if 'transcript_id' in a]
				if len(gs) > 1 or len(tx) > 1:
					print '>1 gene_name or transcript_id found!'
					print gs, tx
					sys.exit(1)
				else:
					gs = gs[0].replace('gene_name', '').replace('\"', '').replace(' ', '')
					tx = tx[0].replace('transcript_id', '').replace('\"', '').replace(' ', '')
					chrom = x[0]
					strand = x[6]
					txstart = int(x[3])
					txend = int(x[4])

					# === get min txstart & max txend per gene symbol ===
					if (gs, chrom, strand) not in gs2minTxStart:
						gs2minTxStart[(gs, chrom, strand)] = txstart
						gs2maxTxEnd[(gs, chrom, strand)] = txend

					elif txend < gs2minTxStart[(gs, chrom, strand)]:  # gene duplication: this transcript is upstream of the previous. keep both
						gs_dup = '_'.join([gs, tx])
						if (gs_dup, chrom, strand) not in gs2minTxStart:
							gs2minTxStart[(gs_dup, chrom, strand)] = txstart
							gs2maxTxEnd[(gs_dup, chrom, strand)] = txend
						elif (gs2minTxStart[(gs_dup, chrom, strand)] != txstart) or (
								gs2maxTxEnd[(gs_dup, chrom, strand)] != txend):
							sys.stderr.write(' '.join([gs_dup, 'exists, but dif transcripts:',
													   str(gs2minTxStart[(gs_dup, chrom, strand)]), str(txstart),
													   str(gs2maxTxEnd[(gs_dup, chrom, strand)]), str(txend)]) + '\n')
							sys.exit(1)

					elif txstart > gs2maxTxEnd[(gs, chrom, strand)]:  # gene duplication: this transcript is downstream of the previous. keep both
						gs_dup = '_'.join([gs, tx])
						if (gs_dup, chrom, strand) not in gs2minTxStart:
							gs2minTxStart[(gs_dup, chrom, strand)] = txstart
							gs2maxTxEnd[(gs_dup, chrom, strand)] = txend
						elif (gs2minTxStart[(gs_dup, chrom, strand)] != txstart) or (
								gs2maxTxEnd[(gs_dup, chrom, strand)] != txend):
							sys.stderr.write(' '.join([gs_dup, 'exists, but dif transcripts:',
													   str(gs2minTxStart[(gs_dup, chrom, strand)]), str(txstart),
													   str(gs2maxTxEnd[(gs_dup, chrom, strand)]), str(txend)]) + '\n')
							sys.exit(1)

					else:  # different isoforms: choose min txstart & max txend
						if txstart < gs2minTxStart[(gs, chrom, strand)]:
							gs2minTxStart[(gs, chrom, strand)] = txstart
						if txend > gs2maxTxEnd[(gs, chrom, strand)]:
							gs2maxTxEnd[(gs, chrom, strand)] = txend

	# write bed file with GS & min txstart & end
	b = open(temp_annotbed, 'w')
	for (gs, chrom, strand) in gs2minTxStart:
		txstart = gs2minTxStart[(gs, chrom, strand)]
		txend = gs2maxTxEnd[(gs, chrom, strand)]
		if ss == 'n':
			b.write('\t'.join(map(str, [chrom, txstart, txend, gs])) + '\n')
		elif ss == 'y':
			b.write('\t'.join(map(str, [chrom, txstart, txend, gs, '0', strand])) + '\n')
	b.close()

	pb.BedTool(temp_annotbed).sort().saveas(temp_annotbed_sort)
	os.remove(temp_annotbed)

	# merge gene annotation
	if ss == 'n':
		pb.BedTool(temp_annotbed_sort).merge(
			c='4', o='distinct', delim=',').saveas(temp_annotbed_merged + '.temp')
	elif ss == 'y':
		pb.BedTool(temp_annotbed_sort).merge(
			c='4,5', o='distinct,distinct', delim=',', s=True).saveas(temp_annotbed_merged + '.temp')

	# reorder fields
	bmo = open(temp_annotbed_merged, 'w')
	with open(temp_annotbed_merged + '.temp', 'r') as f:
		for line in f:
			if ss == 'n':
				(chrom, start, end, name) = line.rstrip().split('\t')
				bmo.write('\t'.join([chrom, start, end, name]) + '\n')
			elif ss == 'y':
				(chrom, start, end, strand, name, score) = line.rstrip().split('\t')
				bmo.write('\t'.join([chrom, start, end, name, score, strand]) + '\n')
	bmo.close()

	check_empty_file(temp_annotbed_merged)
	pb.BedTool(temp_annotbed_merged).sort().saveas(temp_annotbed_merged_sorted)

	os.remove(temp_annotbed_merged + '.temp')
	os.remove(temp_annotbed_merged)
	os.remove(temp_annotbed_sort)


def annotate_tus(infiles, ss, temp_allbed, temp_allbed_sorted, temp_annotbed_merged_sorted,
				 out_intersect, out_merge, out_merge_sort, out_annot):
	"""Label de novo transcription units with gene annotation"""
	# === annotate each TU ===
	o = open(temp_allbed, 'w')
	for infile in infiles:
		f = open(infile, 'r')
		o.write(f.read())
	o.close()

	# sort bed file
	pb.BedTool(temp_allbed).sort().saveas(temp_allbed_sorted)
	os.remove(temp_allbed)

	if ss == 'n':
		pb.BedTool(temp_allbed_sorted).intersect(
			temp_annotbed_merged_sorted, wao=True).saveas(out_intersect)
	elif ss == 'y':
		pb.BedTool(temp_allbed_sorted).intersect(
			temp_annotbed_merged_sorted, wao=True, s=True).saveas(out_intersect)

	os.remove(temp_annotbed_merged_sorted)
	os.remove(temp_allbed_sorted)

	# print '- merging TUs within the same gene'
	gene2minTuStart = {}
	gene2maxTuEnd = {}
	novel_count = 0
	f = open(out_intersect, 'r')
	for line in f:
		x = line.rstrip().split('\t')
		tu_start = int(x[1])
		tu_end = int(x[2])
		if len(x) == 9:  # 4 fields per bed file
			gchrom = x[4]
			if gchrom == '.':
				novel_count += 1
				gene = ':'.join([x[0], x[1], x[2], 'novel' + str(novel_count)])
			else:
				gene = ':'.join(x[4:8])
		elif len(x) == 10:  # 5 fields in first bed file, 4 in 2nd
			gchrom = x[5]
			if gchrom == '.':
				novel_count += 1
				gene = ':'.join([x[0], x[1], x[2], 'novel' + str(novel_count)])
			else:
				gene = ':'.join(x[5:9])
		elif len(x) == 13:  # 6 fields per bed file
			gchrom = x[6]
			if gchrom == '.':
				novel_count += 1
				gene = ':'.join([x[0], x[1], x[2], 'novel' + str(novel_count), x[5]])
			else:
				gene = ':'.join(x[6:10]) + ':' + x[11]
		else:
			sys.stderr.write('EXIT: file type not supported: ' + str(len(x)) + '\n')
			sys.exit(1)

		# merge TUs in the same gene
		if gene not in gene2minTuStart:
			gene2minTuStart[gene] = tu_start
			gene2maxTuEnd[gene] = tu_end
		else:
			if tu_start < gene2minTuStart[gene]:
				gene2minTuStart[gene] = tu_start
			if tu_end > gene2maxTuEnd[gene]:
				gene2maxTuEnd[gene] = tu_end
	f.close()

	# merge any genes in the same TU
	seen = {}
	for gene in gene2minTuStart:
		if len(gene.split(':')) == 4:
			(chrom, start, end, name) = gene.split(':')
			tu = ':'.join(map(str, [chrom, gene2minTuStart[gene], gene2maxTuEnd[gene]]))
		elif len(gene.split(':')) == 5:
			(chrom, start, end, name, strand) = gene.split(':')
			tu = ':'.join(map(str, [chrom, gene2minTuStart[gene], gene2maxTuEnd[gene], strand]))

		if tu not in seen:
			seen[tu] = name
		else:
			seen[tu] = seen[tu] + ',' + name

	# write output
	a = open(out_merge, 'w')
	for tu in seen:
		if len(tu.split(':')) == 3:
			(chrom, start, end) = tu.split(':')
		elif len(tu.split(':')) == 4:
			(chrom, start, end, strand) = tu.split(':')

		name = seen[tu]
		gene_count = 0 if 'novel' in name else len(name.split(','))

		if len(tu.split(':')) == 3:
			a.write('\t'.join(map(str, [chrom, start, end, name, gene_count])) + '\n')
		elif len(tu.split(':')) == 4:
			a.write('\t'.join(map(str, [chrom, start, end, name, gene_count, strand])) + '\n')
	a.close()

	pb.BedTool(out_merge).sort().saveas(out_merge_sort)
	os.remove(out_merge)

	# === merge annotated TUs ===
	if ss == 'n':
		pb.BedTool(out_merge_sort).merge(
			c='4,5', o='distinct,sum', delim=';').saveas(out_annot + '.temp')
	elif ss == 'y':
		pb.BedTool(out_merge_sort).merge(
			c='4,5', o='distinct,sum', delim=';', s=True).saveas(out_annot + '.temp')

	check_empty_file(out_annot + '.temp')

	# reorder fields
	mo = open(out_annot, 'w')
	with open(out_annot + '.temp', 'r') as f:
		for line in f:
			if ss == 'y':
				(chrom, start, end, strand, name, score) = line.rstrip().split('\t')
				mo.write('\t'.join([chrom, start, end, name, score, strand]) + '\n')
			elif ss == 'n':
				(chrom, start, end, name, score) = line.rstrip().split('\t')
				mo.write('\t'.join([chrom, start, end, name, score]) + '\n')

	mo.close()

	os.remove(out_annot + '.temp')
	os.remove(out_merge_sort)
	os.remove(out_intersect)


def get_txregion(in_bed_or_gtf, ss):
	"""Read transcription units into txregions dict"""
	txregions = {}
	f = open(in_bed_or_gtf, 'r')
	for line in f:
		if ss == 'y':
			(chrom, start, end, name, score, strand) = line.rstrip().split('\t')
		elif ss == 'n':
			(chrom, start, end, name, score) = line.rstrip().split('\t')[:5]
			strand = '.'
		start = int(start) + 1  # 0-based -> 1-based
		end = int(end)
		if (chrom, start, end, strand) not in txregions:
			txregions[(chrom, start, end, strand)] = ':'.join([score, name])
		else:
			print 'seen', chrom, start, end, strand
			sys.exit(1)
	f.close()
	return txregions


def get_overlap_edges(tx_start, tx_end, start, end):
	"""Get minimum start and maximum end for each transcrption unit"""
	min_start = 0
	max_end = 0
	if (tx_start >= start) and (tx_start <= end):
		min_start = start
		if (tx_end >= start) and (tx_end <= end):  # 1. txregion is within the gene
			max_end = end
		else:  # 2. txregion overlaps end of gene
			max_end = tx_end
	elif (tx_end >= start) and (tx_end <= end):  # 3. txregion overlaps beginning of gene
		min_start = tx_start
		max_end = end
	elif (start >= tx_start) and (start <= tx_end) and (end >= tx_start) and (end <= tx_end):  # 4. gene is within txregion
		min_start = tx_start
		max_end = tx_end
	return min_start, max_end


def update_gtf_bed(refgtf, ss, txregions, out_annot_gtf, out_annot_bed):
	"""
	Udpate the gtf txstart & txend for each TU that corresponds to an annotated gene.
	Add novel TUs to the reference gtf file.
	Also output a bed file of all TUs.
	"""

	# === gtf file ===
	og = open(out_annot_gtf, 'w')
	seen_txregion = defaultdict(list)
	f = open(refgtf, 'r')
	for line in f:
		if not line.startswith('#'):
			(chrom, source, feature, start, end, score, strand, frame, attr) = line.rstrip().split('\t')
			start = int(start)
			end = int(end)

			if feature == 'gene':
				attr = re.sub(r';$', '', attr)
				attr_list = attr.split('; ')
				attr_gid = [x for x in attr_list if 'gene_id' in x][0]
				attr_gname = [x for x in attr_list if 'gene_name' in x][0]
				gid = attr_gid.replace('gene_id ', '').replace('\"', '')
				gname = attr_gname.replace('gene_name ', '').replace('\"', '')

				# get any txregions that overlap this gene
				overlaps = []
				overlaps.append((start, end))  # include reference start & end, \then get min start & end after getting txregions that overlap
				for (tx_chrom, tx_start, tx_end, tx_strand) in txregions:
					if chrom == tx_chrom:  # get outer edges if they do overlap
						if (ss == 'y' and strand == tx_strand) or ss == 'n':
							min_start, max_end = get_overlap_edges(tx_start, tx_end, start, end)
							if min_start != 0 or max_end != 0:
								overlaps.append((min_start, max_end))
								seen_txregion[(tx_chrom, tx_start, tx_end, tx_strand)].append(gname)

				if len(overlaps) == 0:  # did not overlap any txregions
					og.write(line)
				else:  # get min start & max end if there is more than txregion overlapping this gene
					if len(overlaps) > 1:
						min_start = map(min, zip(*overlaps))[0]
						max_end = map(max, zip(*overlaps))[1]
					else:
						min_start = overlaps[0][0]
						max_end = overlaps[0][1]

					og.write('\t'.join([chrom, source, 'gene', str(min_start), str(max_end),
										score, strand, frame, attr]) + ';\n')

					if min_start != start and max_end != end:
						# add txregion to reference gtf file: transcript
						tx_attr_list = attr_list
						tid = 'transcript_id \"' + '_'.join([chrom, str(min_start), str(max_end), strand, gid]) + '\"'
						tx_attr_list.append(tid)
						tname = 'transcript_name \"' + '_'.join([chrom, str(min_start), str(max_end), strand, gid]) + '\"'
						tx_attr_list.append(tname)
						tx_attr = '; '.join(tx_attr_list)
						og.write('\t'.join([chrom, source, 'transcript', str(min_start), str(max_end),
											score, strand, frame, tx_attr]) + ';\n')

						# add txregion to reference gtf file: exon
						exnum = 'exon_number 1'
						exid = 'exon_id \"' + '_'.join([chrom, str(min_start), str(max_end), strand, gid]) + '\"'
						exon_attr_list = list(tx_attr_list)
						exon_attr_list.extend([exnum, exid])
						exon_attr = '; '.join(exon_attr_list)
						og.write('\t'.join([chrom, source, 'exon', str(min_start), str(max_end),
											score, strand, frame, exon_attr]) + ';\n')
			else:
				og.write(line)
	f.close()

	# === bed file ===
	ob = open(out_annot_bed, 'w')
	novel_count = 1
	for (chrom, start, end, strand) in txregions:
		# get unique gene names only
		(score, name) = txregions[(chrom, start, end, strand)].split(':')
		name_list = []
		for x in name.split(';'):
			for y in x.split(','):
				if y not in name_list:
					name_list.append(y)
		name = ','.join([x for x in name_list if 'novel' not in x])

		# if all annotations novel, just say novel once
		novel_list = [x for x in name_list if 'novel' in x]
		if len(novel_list) == len(name_list):  # all are novel -> rename
			name = 'novel' + str(novel_count)
			novel_count += 1

		# reference gtf: add txregion
		if (chrom, start, end, strand) not in seen_txregion:
			if 'novel' not in name:
				print 'didn\'t overlap gene when it should have?', chrom, start, end, strand, name_list
				sys.exit(1)

			# txregion bed: keep original
			if strand == '.':
				ob.write('\t'.join([chrom, str(start - 1), str(end), name, score]) + '\n')
				gid = 'gene_id \"' + '_'.join([chrom, str(start), str(end)]) + '\"'
				tid = 'transcript_id \"' + '_'.join([chrom, str(start), str(end)]) + '\"'
				tname = 'transcript_name \"' + '_'.join([chrom, str(start), str(end)]) + '\"'

			else:
				ob.write('\t'.join([chrom, str(start - 1), str(end), name, score, strand]) + '\n')
				gid = 'gene_id \"' + '_'.join([chrom, str(start), str(end), str(strand)]) + '\"'
				tid = 'transcript_id \"' + '_'.join([chrom, str(start), str(end), str(strand)]) + '\"'
				tname = 'transcript_name \"' + '_'.join([chrom, str(start), str(end), str(strand)]) + '\"'

			gname = 'gene_name \"' + name + '\"'
			gene_attr = '; '.join([gid, gname])
			tx_attr = '; '.join([gid, tid, gname, tname])

			exnum = 'exon_number 1'
			exid = 'exon_id \"' + '_'.join([chrom, str(end), str(end)]) + '\"'
			exon_attr = '; '.join([gid, tid, gname, tname, exnum, exid])

			if strand == '.':
				og.write('\t'.join([chrom, 'RNA', 'gene', str(start), str(end), '.', '+', '.', gene_attr]) + '\n')
				og.write('\t'.join([chrom, 'RNA', 'transcript', str(start), str(end), '.', '+', '.', tx_attr]) + '\n')
				og.write('\t'.join([chrom, 'RNA', 'exon', str(start), str(end), '.', '+', '.', exon_attr]) + '\n')
			else:
				og.write('\t'.join([chrom, 'RNA', 'gene', str(start), str(end), '.', strand, '.', gene_attr]) + '\n')
				og.write('\t'.join([chrom, 'RNA', 'transcript', str(start), str(end), '.', strand, '.', tx_attr]) + '\n')
				og.write('\t'.join([chrom, 'RNA', 'exon', str(start), str(end), '.', strand, '.', exon_attr]) + '\n')
		elif 'novel' not in name:
			if strand == '.':
				ob.write('\t'.join([chrom, str(start - 1), str(end), name, score]) + '\n')
			else:
				ob.write('\t'.join([chrom, str(start - 1), str(end), name, score, strand]) + '\n')
	og.close()
	ob.close()

	pb.BedTool(out_annot_bed).sort().saveas(out_annot_bed + '.sorted')
	os.rename(out_annot_bed + '.sorted', out_annot_bed)


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									 description='Annotate transcription units and add them to .gtf file')
	parser.add_argument('-i', '--infiles', dest='infiles', type=str, nargs='*', required=True, metavar='',
						help='Input TU bed files from mountainClimberTU output, space-delimited')
	parser.add_argument('-g', '--refgtf', dest='refgtf', type=str, required=True, metavar='',
						help='Reference gtf file')
	parser.add_argument('-s', '--ss', dest='ss', type=str, choices=['y', 'n'], required=True, metavar='',
						help='Strand specific [y, n]')
	parser.add_argument('-o', '--output', dest='output', type=str, required=True, metavar='',
						help='Output prefix')
	args = parser.parse_args()

	# --------------------------------------------------
	# main routine
	# --------------------------------------------------
	print '\njob starting:', str(datetime.now().time())

	# === I/O ===
	temp_allbed = args.output + '.temp.bed'
	temp_allbed_sorted = args.output + '.temp.sorted.bed'

	temp_annotbed = args.output + '.annot.temp'
	temp_annotbed_sort = temp_annotbed + '.sort'
	temp_annotbed_merged = temp_annotbed_sort + '.merged'
	temp_annotbed_merged_sorted = temp_annotbed_sort + '.merged.sort'

	out_intersect = args.output + '.intersect'
	out_merge = args.output + '.preAnnot.bed'
	out_merge_sort = out_merge + '.sort'
	out_annot = args.output + '.annot.bed'

	out_annot_gtf = out_annot.replace('.bed', '') + '.' + os.path.basename(args.refgtf)
	out_annot_bed = out_annot_gtf.replace('.gtf', '.bed')
	out_annot_singleGenes_bed = out_annot_bed.replace('.bed', '_singleGenes.bed')

	# === main ===
	print '- merging overlapping gene annotations', str(datetime.now().time())
	merge_annotation(args.refgtf, args.ss, temp_annotbed, temp_annotbed_sort,
					 temp_annotbed_merged, temp_annotbed_merged_sorted)

	print '- annotating TUs & merging those within the same gene', str(datetime.now().time())
	annotate_tus(args.infiles, args.ss, temp_allbed, temp_allbed_sorted,
				 temp_annotbed_merged_sorted, out_intersect, out_merge, out_merge_sort, out_annot)

	print '- updating gtf to include TUs', str(datetime.now().time())
	if args.ss == 'n':
		print '  -> NOTE: because data is not strand-specific, we assign + strand for all novel TUs in gtf file.'
	txregions = get_txregion(out_annot, args.ss)
	os.remove(out_annot)
	update_gtf_bed(args.refgtf, args.ss, txregions, out_annot_gtf, out_annot_bed)

	print '- getting subset of TUs with <= 1 gene', str(datetime.now().time())
	o = open(out_annot_singleGenes_bed, 'w')
	f = open(out_annot_bed, 'r')
	for line in f:
		score = line.rstrip().split('\t')[4]
		if int(score) <= 1:
			o.write(line)
	f.close()
	o.close()

	print '\nfinished:', str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
