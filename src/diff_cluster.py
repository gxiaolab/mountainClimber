#!/user/bin/python -tt
"""
Cluster change first first across replicates within each condition, and then across conditions.
"""


import os
import sys
import argparse
import numpy as np 					# v1.10.4
import pybedtools as pb
from collections import defaultdict
from datetime import datetime
from sklearn.cluster import DBSCAN 	# v0.18.1


label2priority = {
	'Junction': 1,
	'DistalTSS': 2,
	'DistalPolyA': 3,
	'TandemATSS': 4,
	'TandemAPA': 5,
	'DistalLeft': 6,
	'DistalRight': 7,
	'Exon': 8,
	'Intron': 9,
	'TandemRight': 10,
	'TandemLeft': 11
}


def read_infile(infile, gene2cp, cond, min_fc, min_expn, verbose_flag, lm_flag):
	"""
	Read mountainClimberCP output bed file.
	Get change points with fold change >= min_fc and distal ends.
	For each TU, return:
		a list of change points in the gene meeting the above criteria
		the condition label for this sample
	"""
	count_fc_skipped = 0
	count_total = 0
	with open(infile, 'r') as f:
		for line in f:
			if not line.startswith('track'):
				x = line.rstrip().split('\t')
				chrom = x[0]
				cp_start = x[1]
				cp = x[2]
				name = x[3]
				strand = x[5] if len(x) == 6 else '0'
				gene = ':'.join(name.split(':')[1:4] + [chrom] + [name.split(':')[4]] + [strand])  # gene:gstart:gend:chrom:strand_inferred:strand
				w = int(name.split(':')[5])  # denoise window size

				if lm_flag:
					# clustering diff_cluster output: keep all
					label = name.split(':')[0].split(';')[1]
					gene2cp[gene].append((int(cp), label, w, infile, cond))
				else:
					# clustering mountainClimber output: filter by expression & fold change
					label = name.split(':')[0]
					log2fc = float(x[4])
					fc = round(2 ** abs(log2fc), 2)
					cov_avg = float(name.split(':')[-1])

					if cov_avg >= min_expn and (fc >= min_fc or 'Distal' in label):
						# all change points per gene in all conditions
						gene2cp[gene].append((int(cp), label, w, infile, cond))
					else:
						count_fc_skipped += 1
					count_total += 1

	if verbose_flag:
		print count_fc_skipped, '/', count_total, round(float(count_fc_skipped) / float(count_total) * 100), ' % skipped with fold change <=', min_fc, infile
	return gene2cp


def run_dbscan(w, minpts, cps, labels, infiles, conds):
	"""
	Run DBSCAN on cps with eps=w and min_samples=minpts.
	Return the change points, chang epoint labels, and condition labels for each cluster
	"""
	dbscan = DBSCAN(eps=w, min_samples=minpts).fit_predict(np.asarray(cps).reshape(-1, 1))

	# get change points in each cluster
	db2cps = defaultdict(list)
	db2labels = defaultdict(list)
	db2infiles = defaultdict(list)
	db2conds = defaultdict(list)
	for i, x in enumerate(sorted(cps)):
		if dbscan[i] != -1:
			db2labels[dbscan[i]].append(labels[i])
			db2infiles[dbscan[i]].append(infiles[i])
			db2conds[dbscan[i]].append(conds[i])
			db2cps[dbscan[i]].append(cps[i])

	# get condition label for each cluster: which conditions have >= minpts change points in that cluster
	db2condLabel = {}
	db_delete = []
	for db in db2cps:
		infile_cond_tuple = list(set(zip(db2infiles[db], db2conds[db])))
		condLabel = ''
		for cond in sorted(list(set(db2conds[db]))):
			cond_count = sum(1 for x in infile_cond_tuple if x[1] == cond)
			if cond_count >= minpts:
				condLabel = condLabel + '_' + cond

		if condLabel != '':
			db2condLabel[db] = condLabel[1:]
		else:  # did not reach minpts in any condition
			db_delete.append(db)

	for db in db_delete:
		del db2cps[db]
		del db2labels[db]

	return db2cps, db2labels, db2conds, db2condLabel


def dbscan_each_condition(gene2cp, eps, minpts, ngene2cond_clustered, cond, ss_flag, outfile_cp=0):
	"""
	Run DBSCAN on samples within one condition.
	"""
	if outfile_cp:
		o = open(outfile_cp, 'w')

	total_genes = 0
	for gene in gene2cp:
		if len(gene2cp[gene]) > 0:  # gene may not be expressed in this condition
			cps, labels, ws, infiles, conds = zip(*sorted(gene2cp[gene]))
			if len(cps) == 0:
				continue

			# run dbscan
			w = min(map(int, ws)) / 2 if eps == -1.0 else eps
			db2cps, db2labels, db2conds, db2condLabel = run_dbscan(w, minpts, cps, labels, infiles, conds)

			# write output
			if len(db2cps.keys()) > 0:
				if outfile_cp != 0:  # clustering diff_cluster output -> do not output
					cp_med_list, cp_med2label = write_cp(gene, w, db2cps, db2labels, db2condLabel, o, ss_flag)
				total_genes += 1
				ngene2cond_clustered[gene].append(cond)

	if outfile_cp != 0:
		o.close()

	return ngene2cond_clustered, total_genes


def dbscan_across_conditions(cond2cpfile, cond2gene2cp_all, ngene2cond_clustered, min_conditions,
	eps, cond2minpts, ss_flag, outfile_cp, outfile_seg):
	"""
	Run DBSCAN on each gene across conditions.
	"""
	cond2gene_clustered = {}
	o = open(outfile_cp, 'w')
	t = open(outfile_seg, 'w')

	all_genes = set(list([x for y in cond2cpfile.keys() for x in sorted(cond2gene2cp_all[y].keys())]))
	for gene in all_genes:
		if len(ngene2cond_clustered[gene]) >= min_conditions:
			all_cps = [x for y in cond2cpfile.keys() for x in sorted(cond2gene2cp_all[y][gene])]
			if len(all_cps) == 0:
				continue

			cps, labels, ws, infiles, conds = zip(*all_cps)

			# run dbscan
			w = min(map(int, ws)) / 2 if eps == -1.0 else eps
			db2cps, db2labels, db2conds, db2condLabel = run_dbscan(w, min(cond2minpts.values()),
				cps, labels, infiles, conds)

			# filter change points that did not meet args.minpts cutoff in min_conditions total conditions
			db_delete = []
			for d, db in enumerate(db2cps.keys()):
				conds_unique = list(set(db2conds[db]))
				if len(conds_unique) == 1:
					cond = conds_unique[0]
					if db2conds[db].count(cond) < cond2minpts[cond]:
						db_delete.append(db)
				else:
					db_keep_cond_count = 0
					for cond in conds_unique:
						if db2conds[db].count(cond) >= cond2minpts[cond]:
							db_keep_cond_count += 1
					if db_keep_cond_count < min_conditions:
						db_delete.append(db)

			for db in list(set(db_delete)):
				del db2cps[db]
				del db2labels[db]

			# write output
			if len(db2cps.keys()) > 0:
				cp_med_list, cp_med2label = write_cp(gene, w, db2cps, db2labels, db2condLabel, o, ss_flag)
				write_seg(gene, w, cp_med_list, cp_med2label, t, ss_flag)

				for cond in ngene2cond_clustered[gene]:
					if cond not in cond2gene_clustered:
						cond2gene_clustered[cond] = 1
					else:
						cond2gene_clustered[cond] += 1
	o.close()
	t.close()
	return cond2gene_clustered


def write_cp(gene, w, db2cps, db2label, db2condLabel, o, ss_flag):
	"""
	Write output: bed file of change points.
	"""
	cp_med_list = []
	cp_med2label = {}

	(gs, gstart, gend, chrom, inferred_strand, strand) = gene.split(':')
	out_gene = ':'.join(gene.split(':')[:-1])  # exclude "0" if non-strand-specific
	dbmax = max(db2cps.keys())
	for d, db in enumerate(db2cps):
		cp_med = int(np.median(db2cps[db]))
		cp_med_list.append(cp_med)
		cp_std = round(np.std(db2cps[db]), 2)
		labels = sorted(list(set(db2label[db])))
		label = ','.join(sorted(list(set(db2label[db]))))
		condLabel = db2condLabel[db]

		label_priorities = [label2priority[x] for x in labels]
		label_prioritized = labels[label_priorities.index(min(label_priorities))]

		if cp_med not in cp_med2label:
			cp_med2label[cp_med] = (label, label_prioritized, condLabel, db, dbmax)
		else:
			sys.stderr.write('EXIT: two clusters have the same median value!\n')
			sys.stderr.write('\t'.join(map(str, [cp_med, label, db, dbmax, cp_med2label[cp_med]])) + '\n')
			sys.exit(1)

		# all change points
		if ss_flag:
			o.write('\t'.join(map(str, [chrom, cp_med - 1, cp_med, ':'.join(map(str,
				[label_prioritized + ';' + condLabel, out_gene, w, min(db2cps[db]),
				max(db2cps[db]), cp_std, len(db2cps[db])])), d, strand])) + '\n')
		else:
			o.write('\t'.join(map(str, [chrom, cp_med - 1, cp_med, ':'.join(map(str,
				[label_prioritized + ';' + condLabel, out_gene, w, min(db2cps[db]),
				max(db2cps[db]), cp_std, len(db2cps[db])])), d])) + '\n')

	return cp_med_list, cp_med2label


def write_seg(gene, w, cp_med_list, cp_med2label, t, ss_flag):
	"""
	Write output: bed file of segments (before & after each change point)
	"""
	(gs, gstart, gend, chrom, inferred_strand, strand) = gene.split(':')
	out_gene = ':'.join(gene.split(':')[:-1])  # exclude "0" if non-strand-specific
	cp_med_list_sorted = sorted(cp_med_list)
	for i, cp in enumerate(cp_med_list_sorted):
		if i > 0:
			(label, label_prioritized, condLabel, db, dbmax) = cp_med2label[cp]
			(prev_label, prev_label_prioritized,
				prev_condLabel, prev_db, prev_dbmax) = cp_med2label[cp_med_list_sorted[i - 1]]

			if ss_flag:
				t.write('\t'.join(map(str, [chrom, cp_med_list_sorted[i - 1], cp,
					';'.join([prev_label_prioritized, prev_condLabel]) + '|' + ';'.join(
						[label_prioritized, condLabel]) + ':' + out_gene + ':' + str(w), i, strand])) + '\n')
			else:
				t.write('\t'.join(map(str, [chrom, cp_med_list_sorted[i - 1], cp,
					';'.join([prev_label_prioritized, prev_condLabel]) + '|' + ';'.join(
						[label_prioritized, condLabel]) + ':' + out_gene + ':' + str(w), i])) + '\n')


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Cluster change points first across replicates within each condition, and then across conditions.')
	group = parser.add_argument_group('Input')
	group.add_argument('-i', '--input', dest='input', type=str, nargs='*', metavar = '',
		help='List of space-delimited change point files.')
	group.add_argument('-c', '--conditions', dest='conditions', type=str, nargs='*', metavar = '',
		help='List of space-delimited condition labels for each --input file.')

	group = parser.add_argument_group('Parameters')
	group.add_argument('-n', '--minpts', dest='minpts', type=int, nargs='*', metavar = '',
		help='List of space-delimited DBSCAN minPts values. These indicate the minimum # points for \
		DBSCAN to consider a core point. The minimum of this list will be used to cluster across conditions.')
	group.add_argument('-e', '--eps', dest='eps', type=float, default=-1.0, metavar = '',
		help='Maximum distance between 2 points in a neighborhood. \
		-1.0 indicates using the minimum optimal window size from mountain climber.')
	group.add_argument('-f', '--min_fc', dest='min_fc', type=float, default=-1.0, metavar = '',
		help='Minimum fold change for change points.')
	group.add_argument('-d', '--min_conditions', dest='min_conditions', type=float, default=1, metavar = '',
		help='Minimum number of conditions for a gene to be clustered across conditions.')
	group.add_argument('-m', '--min_expn', dest='min_expn', type=float, default=0, metavar = '',
		help='Minimum expression in exons for a gene to be clustered.')
	group.add_argument('-l', '--lm_flag', dest='lm_flag', action='store_true',
		help='Input are results from diff_cluster.')
	group.add_argument('-s', '--ss_flag', dest='ss_flag', action='store_true',
		help='Flag: RNA-Seq is strand-specific.')

	group = parser.add_argument_group('output')
	group.add_argument('-o', '--output', dest='output', type=str, metavar = '',
		help='Output prefix. Output files include _cluster_totals.txt, _segments.bed, _cp.bed, and one _cp.bed file for each condition. _cp.bed name field = label_prioritized;condition_labels:gene:TUstart:TUend:chrom:strand:dbscan_epsilon:min_clustered_change_point:max_clustered_change_point:cluster_standard_deviation:total_clusters. _segments.bed name field = label_prioritized_cp1;condition_labels_cp1|label_prioritized_cp2;condition_labels_cp2:gene:TUstart:TUend:chrom:strand:dbscan_epsilon.')
	group.add_argument('-v', '--verbose', dest='verbose', action='store_true',
		help='Print progress details.')
	args = parser.parse_args()

	if args.verbose:
		print args

	# --------------------------------------------------
	# main routine
	# --------------------------------------------------
	print '\njob starting:', str(datetime.now().time())

	if not args.input:
		sys.stderr.write('EXIT: Please provide --input')
		sys.exit(1)

	if not args.output:
		sys.stderr.write('EXIT: Please provide --output')
		sys.exit(1)

	if not args.minpts:
		sys.stderr.write('WARNING: no minpts provided: setting equal to 1\n')
		args.minpts = [1] * len(args.conditions)

	if len(args.input) != len(args.conditions) or len(args.input) != len(args.minpts):
		sys.stderr.write(' '.join(map(str, ['EXIT: number of samples don\'t match!:',
			len(args.input), 'input,', len(args.conditions), 'conditions,', len(args.minpts), 'minpts'])) + '\n')
		sys.exit(1)

	# === I/O ===
	outfile_cp = args.output + '_cp.bed'
	outfile_seg = args.output + '_segments.bed'
	outfile_summary = args.output + '_cluster_totals.txt'

	# === get the input files for each condition ===
	cond2cpfile = defaultdict(list)
	cond2minpts = {}
	for c, cond in enumerate(args.conditions):
		cond2cpfile[cond].append(args.input[c])
		if cond not in cond2minpts:
			cond2minpts[cond] = args.minpts[c]
		elif cond2minpts[cond] != args.minpts[c]:
			sys.stderr.write(' '.join(map(str, ['EXIT: 2 different minpts provided for', cond, ':', args.minpts[c], cond2minpts[cond]])) + '\n')
			sys.exit(1)

	# === (1) filter change points by fold change, (2) filter TUs by reproducibility ===
	print '- reading input change point files', str(datetime.now().time())
	cond2gene2cp_all = {}
	cond2nGenes = {}
	cond2nGenesFiltered = {}
	for cond in sorted(cond2cpfile.keys()):
		print '  -', cond, str(datetime.now().time())

		# === get all change points with >= min_fc ===
		gene2cp = defaultdict(list)
		for cpfile in cond2cpfile[cond]:
			gene2cp = read_infile(cpfile, gene2cp, cond, args.min_fc,
				args.min_expn, args.verbose, args.lm_flag)
		cond2nGenes[cond] = len(gene2cp.keys())  # count total genes
		if args.verbose:
			print ' -> total genes:', len(gene2cp.keys())

		# === all genes across all conditions ===
		cond2gene2cp_all[cond] = defaultdict(list)
		for gene in gene2cp:
			cond2gene2cp_all[cond][gene].extend(gene2cp[gene])

	# === DBSCAN ===
	print '- DBSCAN: in each condition', str(datetime.now().time())
	cond2ngenes_eachCondition = {}
	ngene2cond_clustered = defaultdict(list)
	for cond in sorted(cond2cpfile.keys()):
		print '  -', cond, str(datetime.now().time())
		if len(cond2cpfile.keys()) == 1:  # only one condition -> print output bed file
			ngene2cond_clustered, total_genes = dbscan_each_condition(cond2gene2cp_all[cond],
				args.eps, cond2minpts[cond], ngene2cond_clustered, cond, args.ss_flag,
				outfile_cp=outfile_cp)
		else:  # >1 conditions -> write condition-specific outputs separately
			outfile_cp_cond = args.output + '_cp_' + cond + '.bed'
			ngene2cond_clustered, total_genes = dbscan_each_condition(cond2gene2cp_all[cond],
				args.eps, cond2minpts[cond], ngene2cond_clustered, cond, args.ss_flag,
				outfile_cp=outfile_cp_cond)

		cond2ngenes_eachCondition[cond] = total_genes

	# === dbscan across conditions: keep outliers ===
	print '- DBSCAN: across conditions', str(datetime.now().time())
	if len(cond2cpfile.keys()) > 1:
		cond2gene_clustered = dbscan_across_conditions(cond2cpfile, cond2gene2cp_all,
			ngene2cond_clustered, args.min_conditions, args.eps, cond2minpts,
			args.ss_flag, outfile_cp=outfile_cp, outfile_seg=outfile_seg)

		# === print totals across conditions ===
		o = open(outfile_summary, 'w')
		o.write('\t'.join(['Condition', 'Total Samples', 'DBSCAN minPts', 'Total Genes',
			'Total genes clustered per condition', 'Total genes clustered across conditions']) + '\n')
		for cond in cond2nGenes:
			if cond not in cond2gene_clustered:
				cond2gene_clustered[cond] = 'NA'
			o.write('\t'.join(map(str, [cond, len(cond2cpfile[cond]), cond2minpts[cond], cond2nGenes[cond], cond2ngenes_eachCondition[cond], cond2gene_clustered[cond]])) + '\n')
		o.close()
	else:
		print '  -> only one condition found'

		# === print totals across conditions ===
		o = open(outfile_summary, 'w')
		o.write('\t'.join(['Condition', 'Total Samples', 'DBSCAN minPts', 'Total Genes',
			'Total genes clustered per condition']) + '\n')
		for cond in cond2nGenes:
			o.write('\t'.join(map(str, [cond, len(cond2cpfile[cond]), cond2minpts[cond], cond2nGenes[cond], cond2ngenes_eachCondition[cond]])) + '\n')
		o.close()

	# === sort output ===
	pb.BedTool(outfile_cp).sort().saveas(outfile_cp + '.sorted')
	os.rename(outfile_cp + '.sorted', outfile_cp)
	print '\nfinished:', str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
