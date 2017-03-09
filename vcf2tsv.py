#!/usr/bin/env python2
"""
    Transform a VCF file into a tab-delimited file with 1s and 0s for
    genotype calls.

    vcf2tsv.py --output <output file> <input file>

    OPTIONS

        --subset        A bed file with positions to include, default is
                        to include everything
        --rs            Include SNP ID
        --min-gt        Minimum genotyping rate (0)
        --impute        Replace missing genotypes with the mean across
                        the entire sample
"""

#==============================================================================#

import argparse
import re

import vcf

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--subset')
parser.add_argument('--min-gt', default=0.0, type=float)
parser.add_argument('--impute', default=False, action='store_true')
parser.add_argument('--rs', default=False, action='store_true')
parser.add_argument('input')
args = parser.parse_args()

if args.subset is None:
    positions = None
else:
    positions = set()
    with open(args.subset, 'rb') as handle:
        for line in handle:
            contig, _, pos = line.strip().split()[:3]
            positions.add((contig, pos))

with open(args.output, 'wb') as out:
    rdr = vcf.Reader(filename=args.input)
    out.write('contig\tpos\t' + '\t'.join(rdr.samples) + '\n')
    for rec in rdr:
        if (positions is not None) and ((rec.CHROM, str(rec.POS)) not in positions):
            continue
        rs = rec.ID
        gts = []
        missed = 0
        sum_gt = 0
        for s in rec.samples:
            gt = re.split('\||/', s['GT'])
            if gt[0] == '.':
                g = 'NA'
                missed += 1
            else:
                g = sum(1 for c in gt if c == '0') / 2.
                sum_gt += g
            gts.append(str(g))
        if 1 - (missed / float(len(gts))) < args.min_gt:
            continue
        if missed > 0 and args.impute:
            mean_gt = sum_gt / float(len(gts) - missed)
            gts = [g if g != 'NA' else mean_gt for g in gts]
        out.write(rec.CHROM + '\t' + str(rec.POS))
        if args.rs:
            out.write('\t' + rs)
        out.write('\t' + '\t'.join(str(g) for g in gts) + '\n')
