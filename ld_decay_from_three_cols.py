#!/usr/bin/env python2.7
"""
    Calculate LD decay given a file with variant IDs and R^2 values.

    ld_decay_from_three_cols.py <ld file> <plink bim file> <replicon> <replicon length>
"""

import argparse
import gzip
import math
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--limit', type=int, default=500000)
parser.add_argument('ld')
parser.add_argument('bim')
parser.add_argument('replicon')
parser.add_argument('length', type=int)
args = parser.parse_args()

target = args.replicon
target_len = args.length
max_dist = args.limit
half_replicon = math.ceil(target_len / 2.)

positions = {}
with open(args.bim, 'rb') as ih:
    for line in ih:
        chrom, rs, _, pos, _, _ = line.strip().split()
        positions[rs] = (chrom, int(pos))

ofun = open
if args.ld.endswith('.gz'):
    ofun = gzip.open
with ofun(args.ld, 'rb') as ih:
    ih.readline()
    for line in ih:
        rs0, rs1, r2 = line.strip().split('\t')
        p0 = positions[rs0]
        p1 = positions[rs1]
        if p0[0] != target or p1[0] != target:
            continue
        d = abs(p1[1] - p0[1])
        if d > half_replicon:
            d = target_len - d
        if d > max_dist:
            continue
        sys.stdout.write(r2 + '\t' + str(d) + '\t' + str(math.ceil(d / 1000.)) + '\n')
