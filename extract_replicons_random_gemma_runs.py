#!/usr/bin/env python3
"""
    extract_replicons_random_gemma_runs.py --output <output file>
        [-n (50)]
        <annotated file> <combined random file>
"""

import argparse
import collections
import csv
import gzip
import itertools

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('-n', default=50, type=int)
parser.add_argument('annot')
parser.add_argument('combined')
args = parser.parse_args()

n = args.n

op = open
if args.annot.endswith('.gz'):
    op = gzip.open
replicons = collections.defaultdict(list)
with op(args.annot, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        replicons[row[3]].append(row[0])

op = open
if args.combined.endswith('.gz'):
    op = gzip.open

with open(args.output, 'wt') as oh:
    oh.write('run\trank\tgroup\treplicon\n')
    with op(args.combined, 'rt') as ih:
        rdr = csv.DictReader(ih, delimiter='\t')
        for run, rows in itertools.groupby(rdr, lambda x: x['run']):
            ldgroups = []
            for row in rows:
                ldgroups.append((float(row['p_lrt']), row['rs']))
            ldgroups.sort(key=lambda x: x[0])
            for i in range(n):
                for r in replicons[ldgroups[i][1]]:
                    oh.write(str(run) + '\t' + str(i) + '\t' + ldgroups[i][1] + '\t'
                             + r + '\n')
