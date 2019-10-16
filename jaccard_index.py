#!/usr/bin/env python3
"""
    Calculate the Jaccard index given a tsv file with a column containing
    comma-separated lists of items.

    jaccard_index.py <input file> <column>
"""

import argparse
import collections
import csv
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('input')
parser.add_argument('column')
args = parser.parse_args()

shared = collections.defaultdict(int)
total = collections.defaultdict(int)
with open(args.input, 'rt') as ih:
    rdr = csv.DictReader(filter(lambda x: not x.startswith('#'), ih), delimiter='\t')
    for i, row in enumerate(rdr):
        groups = row[args.column].split(',')
        for g1 in groups:
            total[g1] += 1
            for g2 in groups:
                shared[(g1, g2)] += 1

names = sorted(total.keys())
sys.stdout.write('\t' + '\t'.join(names) + '\n')
for n in names:
    sys.stdout.write(n)
    for n2 in names:
        j = shared[(n, n2)] / float(total[n] + total[n2] - shared[(n, n2)])
        sys.stdout.write('\t' + str(j))
    sys.stdout.write('\n')
