#!/usr/bin/env python3
"""
    Given GEMMA LMM output (with LRT test) from just on representative of
    each LD group, add in all the other variants associated with that
    LD group.
"""

import argparse
import csv
import os
import os.path as osp
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output-bed')
parser.add_argument('gemma')
parser.add_argument('onevar')
args = parser.parse_args()

ld_groups = dict()
with open(args.onevar, 'r') as ih:
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

rows = []
with open(args.gemma, 'r') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    colnames = rdr.fieldnames
    for row in rdr:
        rows.append(row)
rows.sort(key=lambda x: float(x['p_lrt']))

with open(args.output_bed, 'w') as oh:
    for i, row in enumerate(rows):
        for var in ld_groups[row['rs']]:
            oh.write('\t'.join(str(x) for x in var[:3]) + '\t' +
                     row['rs'] + '\t' + str(i) + '\n')
