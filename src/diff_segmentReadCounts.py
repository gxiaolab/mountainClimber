#!/user/bin/python -tt
"""
Calculate the average reads/bp for each segment after clustering.
"""


import os
import sys
import argparse
import numpy as np  # v1.10.4
import pybedtools as pb
from collections import defaultdict
from datetime import datetime


def bedgraph_per_gene_ss(genes, bg_plus, bg_minus, bgfile):
    """
    Run bedtools intersect: strand-specific.
    Keep strands separate so that lines of coverage for a given gene are consecutive.
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

    # === bedtools intersect ===
    pb.BedTool(plus_bed).sort().intersect(bg_plus, wo=True, sorted=True).saveas(bgfile + '.plus')
    pb.BedTool(minus_bed).sort().intersect(bg_minus, wo=True, sorted=True).saveas(bgfile + '.minus')

    t = open(bgfile, 'w')
    t.write(open(bgfile + '.plus').read())
    t.write(open(bgfile + '.minus').read())
    t.close()

    for file in [bgfile + '.plus', bgfile + '.minus', plus_bed, minus_bed]:
        os.remove(file)


def bedgraph_per_gene_nss(genes, bg, bgfile):
    """Run bedtools intersect: non-strand-specific"""
    pb.BedTool(genes).sort().intersect(bg, wo=True, sorted=True).saveas(bgfile)


def get_seg2cov(intersect, cond, sample, outfile):
    """Get coverage of each segment & write to outfile in bed format"""
    seg2cov = {}
    o = open(outfile, 'w')

    maxl = 0
    with open(intersect, 'r') as f:
        for l, line in enumerate(f):
            maxl = l

    with open(intersect, 'r') as f:
        for l, line in enumerate(f):
            if line != '':
                x = line.rstrip().split('\t')
                if len(x) == 11:
                    (achrom, astart, aend, ageneid, ascore, astrand, bchrom, bstart, bend, bcov, overlap_len) = x
                elif len(x) == 10:
                    (achrom, astart, aend, ageneid, ascore, bchrom, bstart, bend, bcov, overlap_len) = x
                    astrand = 0
                else:
                    sys.stderr.write('EXIT: do not recognize bedgraph intersect format\n')
                    sys.stderr.write(line)
                    sys.exit(1)

                astart = int(astart)
                aend = int(aend)
                bstart = int(bstart)
                bend = int(bend)
                bcov = float(bcov)
            else:
                x = ''

            if overlap_len == '0':
                continue

            if l == 0:  # first line
                prev_gene = ':'.join(x[:5]) if astrand == 0 else ':'.join(x[:6])
                this_start = max(astart, bstart)
                this_end = min(aend, bend)
                prev_cov_array = np.zeros(aend - astart)
                prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

                # === next round  ===
                prev_start = astart
                prev_end = aend
                prev_geneid = ageneid
                prev_chrom = achrom
                prev_strand = astrand

                if l == maxl:
                    cov_avg = np.mean(prev_cov_array)
                    if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
                        seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
                        if prev_strand == 0:
                            o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
                        else:
                            o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
                    else:
                        sys.stderr.write(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
                        sys.stderr.write(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
                        sys.exit(1)
            else:
                this_gene = ':'.join(x[:5]) if astrand == 0 else ':'.join(x[:6])
                if line == '' and this_gene == prev_gene and this_gene == '':  # EOF
                    break
                elif this_gene == prev_gene and l != maxl:  # get coverage
                    this_start = max(astart, bstart)
                    this_end = min(aend, bend)
                    prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

                    # === next round  ===
                    prev_start = astart
                    prev_end = aend
                    prev_geneid = ageneid
                    prev_chrom = achrom
                    prev_strand = astrand
                else:  # finished reading all info for one gene
                    # === get per-base coverage ===
                    if l == maxl:
                        # === last line of this gene & last line of the file ===
                        this_start = max(astart, bstart)
                        this_end = min(aend, bend)

                        # === previous gene ===
                        prev_cov_array[(this_start - astart):(this_end - astart)] += bcov
                        cov_avg = np.mean(prev_cov_array)
                        if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
                            seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
                            if prev_strand == 0:
                                o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
                            else:
                                o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
                        else:
                            sys.stderr.write(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
                            sys.stderr.write(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
                            sys.exit(1)

                        if this_gene != prev_gene:
                            prev_start = astart
                            prev_end = aend
                            prev_geneid = ageneid
                            prev_chrom = achrom
                            prev_strand = astrand

                            cov_avg = np.mean(prev_cov_array)
                            if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
                                seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
                                if prev_strand == 0:
                                    o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
                                else:
                                    o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
                            else:
                                sys.stderr.write(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
                                sys.stderr.write(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
                                sys.exit(1)

                    else:
                        cov_avg = np.mean(prev_cov_array)
                        if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
                            seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
                            if prev_strand == 0:
                                o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
                            else:
                                o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
                        else:
                            sys.stderr.write(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
                            sys.stderr.write(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
                            sys.exit(1)

                        # === first line of the next gene ===
                        this_start = max(astart, bstart)
                        this_end = min(aend, bend)
                        prev_cov_array = np.zeros(aend - astart)
                        prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

                        # === next round ===
                        prev_start = astart
                        prev_end = aend
                        prev_geneid = ageneid
                        prev_chrom = achrom
                        prev_strand = astrand
                        prev_gene = this_gene
    o.close()


def main(argv):
    # --------------------------------------------------
    # get args
    # --------------------------------------------------
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Calculate the average reads/bp for each segment after clustering.')
    group = parser.add_argument_group('Input')
    group.add_argument('-i', '--segments', dest='segments', type=str, metavar='', help='_segments.bed output from clusterCPs')
    group.add_argument('-p', '--bgplus', dest='bgplus', type=str, nargs='*', metavar='', help='List of space-delimited bedgraphs: non-strand-specific or plus strand.')
    group.add_argument('-m', '--bgminus', dest='bgminus', type=str, nargs='*', metavar='', help='List of space-delimited bedgraphs: minus strand.')
    group.add_argument('-c', '--conditions', dest='conditions', type=str, nargs='*', metavar='', help='List of space-delimited condition labels for each --bgplus file')
    group = parser.add_argument_group('Output')
    group.add_argument('-o', '--output', dest='output', type=str, metavar='', help='Output prefix.')
    args = parser.parse_args()
    print args

    # --------------------------------------------------
    # main routine
    # --------------------------------------------------
    print '\njob starting:', str(datetime.now().time())

    if not args.segments:
        sys.stderr.write('EXIT: Please provide --segments')
        sys.exit(1)

    if not args.conditions:
        sys.stderr.write('EXIT: Please provide --conditions')
        sys.exit(1)

    if args.bgplus:
        if len(args.conditions) != len(args.bgplus):
            sys.stderr.write('\n'.join(['EXIT: number of samples don\'t match!:'] + '\n'))
    else:
        sys.stderr.write('EXIT: Please provide --bgplus')
        sys.exit(1)

    if args.bgminus:
        if len(args.conditions) != len(args.bgminus):
            sys.stderr.write('\n'.join(['EXIT: number of samples don\'t match!:'] + '\n'))

    # === set temporary dir ===
    tmpdir = args.output
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    pb.set_tempdir(tmpdir)

    # === get the input files for each condition ===
    cond2bgplus = defaultdict(list)
    cond2bgminus = defaultdict(list)
    for c, cond in enumerate(args.conditions):
        cond2bgplus[cond].append(args.bgplus[c])
        if args.bgminus:
            cond2bgminus[cond].append(args.bgminus[c])

    # === get bedgraph for each segment ===
    print '- getting read counts per bp for each gene', str(datetime.now().time())
    # condSeg2cov = {}
    for cond in cond2bgplus:
        for b, bgplus in enumerate(cond2bgplus[cond]):
            sample = os.path.splitext(os.path.basename(bgplus))[0]
            print '  -> bedtools intersect:', cond, sample, str(datetime.now().time())

            # get bedgraph for each segment
            intersect = '_'.join([args.output, cond, sample, 'intersect.txt'])
            if args.bgminus:
                bedgraph_per_gene_ss(args.segments, bgplus, cond2bgminus[cond][b], intersect)
            else:
                bedgraph_per_gene_nss(args.segments, bgplus, intersect)

            print '  - avg coverage per bp for each segment', str(datetime.now().time())
            get_seg2cov(intersect, cond, sample, outfile='_'.join([args.output, sample, cond + '_readCounts.bed']))

            # remove temporary file
            os.remove(intersect)

    # if os.listdir(tmpdir) == []:
    #   os.rmdir(tmpdir)
    # else:
    #   sys.stderr.write('WARNING: could not remove directory ' + tmpdir + '\n')

    print '\nfinished:', str(datetime.now().time())


# boilerplate
if __name__ == '__main__':
    main(sys.argv[1:])
