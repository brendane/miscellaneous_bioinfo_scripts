#!/usr/bin/env python2.7
"""
    Take a gemma -lmm 4 output file that was created with LD groups and
    make a bed file with individual variants.

    This is a modified version of merge_gemma_lmm_iters.py.

    Note that if the position is > 5000000, it will subtract 5 Mb; this
    is the way I coded single-reference CNVs.
"""

P_COLUMN = 'p_lrt'

import argparse
import csv
import gzip
import os
import os.path as osp
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('infile')
parser.add_argument('onevar')
args = parser.parse_args()

# Read the LD grouping information
ld_groups = dict()
with open(args.onevar, 'rb') as ih:
    for line in ih:
        row = line.strip().split('\t')
        group = 'group-' + row[4]
        vars_in_group = row[5].split(',')
        vars_list = []
        for v in vars_in_group:
            try:
                t, _, c, p, _ = v.split('-')
            except:
                t, c, p, _ = v.split('-')
            p = int(p)
            if p > 5000000:
                p -= 5000000
            vars_list.append((c, p-1, p, v))
        ld_groups[group] = vars_list

assoc_file = args.infile

# Read the association results and sort by p-value.
rows = []
ofun = open
if args.infile.endswith('.gz'):
    ofun = gzip.open
with ofun(args.infile, 'rb') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        rows.append(row)
rows.sort(key=lambda x: float(x[P_COLUMN]))

# Write the output file
with open(args.output, 'wb') as oh:
    for i, row in enumerate(rows):
        for var in ld_groups[row['rs']]:
            oh.write('\t'.join(str(x) for x in var[:3]) + '\t' +
                     row['rs'] + '\t' + str(i) + '\n')
