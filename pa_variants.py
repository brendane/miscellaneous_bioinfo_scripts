#!/usr/bin/env python2
"""
    Tabulate presence absence variants in plink tped format.

    pa_variants.py [options] --output <output prefix> <input file>

        --lower         Less than or equal to lower propotion of sites
                        covered means a feature is absent (0.5)
        --upper         Greater than or equal to upper means present (0.5)
        --pad           Add bp to the position so that there won't be any
                        overlap with SNP variants (0)
"""

import argparse
import csv
import os
import os.path as osp
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--lower', type=float, default=0.5)
parser.add_argument('--upper', type=float, default=0.5)
parser.add_argument('--pad', type=int, default=0)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

lower_cov = args.lower
upper_cov = args.upper
pad = args.pad
out_prefix = args.output

covered = {}
with open(args.input, 'rb') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    strains = sorted(s for s in rdr.fieldnames if s not in {'feature', 'replicon', 'pos'})
    for row in rdr:
        feature = row['feature']
        replicon = row['replicon']
        pos = str(int(row['pos']) + pad)
        covd = {}
        for strain in strains:
            c = float(row[strain])
            if c <= lower_cov:
                covd[strain] = 'A'
            elif c >= upper_cov:
                covd[strain] = 'T'
            else:
                covd[strain] = '0'
        covered[(feature, replicon, pos)] = covd

strains.sort()
with open(out_prefix + '.tfam', 'wb') as out:
    for strain in strains:
        out.write(strain + '\t' + strain + '\t0\t0\t0\t-9\n')

with open(out_prefix + '.tped', 'wb') as out:
    for feature in covered:
        out.write(feature[1] + '\t' + feature[0] + '-' + feature[2] + 
                  '\t0\t' + feature[2])
        for strain in strains:
            gt = covered[feature][strain]
            out.write('\t' + gt + '\t' + gt)
        out.write('\n')
