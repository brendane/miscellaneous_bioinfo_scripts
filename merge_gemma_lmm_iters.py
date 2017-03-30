#!/usr/bin/env python2.7
"""
    Merge a set of LMM runs from GEMMA. Assumes the runs were done in
    the same way that the 2017-03-19 metab_gwas runs were done.
"""

import argparse
import csv
import gzip
import os
import os.path as osp
import re

def get_top_var(rdr, itern):
    best_p = 1.
    best_row = None
    for row in rdr:
        p = float(row['p_lrt'])
        if p < best_p:
            best_p = p
            best_row = row
    best_row['iter'] = itern
    return best_row

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--output-bed')
parser.add_argument('indir')
parser.add_argument('onevar')
args = parser.parse_args()

ld_groups = dict()
with open(args.onevar, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
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

assoc_files = map(lambda x: (int(re.sub('\\.assoc.+', '', re.sub('output', '', x))), x),
                  filter(lambda x: x.endswith('assoc.txt.gz'),
                         os.listdir(args.indir)))
assoc_files.sort(key=lambda x: x[0])
last_run = max(x[0] for x in assoc_files)

top_vars = []
for i in range(last_run):
    with gzip.open(osp.join(args.indir, assoc_files[i][1]), 'rb') as ih:
        rdr = csv.DictReader(ih, delimiter='\t')
        top_vars.append((i, get_top_var(rdr, i)))

rows = []
with gzip.open(args.output, 'wb') as oh:
    with gzip.open(osp.join(args.indir, assoc_files[last_run][1]), 'rb') as ih:
        rdr = csv.DictReader(ih, delimiter='\t')
        colnames = ['iter'] + rdr.fieldnames
        oh.write('\t'.join(colnames) + '\n')
        for i, row in top_vars:
            rows.append(row)
            oh.write(str(i))
            for cn in colnames[1:]:
                oh.write('\t' + row[cn])
            oh.write('\n')
        for row in rdr:
            row['iter'] = last_run
            rows.append(row)
    rows.sort(key=lambda x: (int(x['iter']), float(x['p_lrt'])))
    for row in rows:
        for c, cn in enumerate(colnames):
            if c > 0:
                oh.write('\t' + str(row[cn]))
            else:
                oh.write(str(row[cn]))
        oh.write('\n')

with open(args.output_bed, 'wb') as oh:
    for i, row in enumerate(rows):
        for var in ld_groups[row['rs']]:
            oh.write('\t'.join(str(x) for x in var[:3]) + '\t' +
                     row['rs'] + '\t' + str(i) + '\n')
