#!/usr/bin/env python3
"""
    1. Filter out 1-M and M-M orthologs
    2. Filter out orthologs that differ by more than 5% in length
    3. Extract locus tags
"""

import argparse
import collections
import csv
import itertools
import sys

import networkx

parser = argparse.ArgumentParser()
parser.add_argument('--threshold', type=float)
parser.add_argument('--use-third')
parser.add_argument('files', nargs='+')
args = parser.parse_args()

thresh = args.threshold
files = args.files
if args.use_third is None:
    third = []
else:
    third = args.use_third.split(',')

g = networkx.Graph()
strains = set()
for fname in files:
    with open(fname, 'r') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        s0, s1 = next(rdr)[2:4]
        for row in rdr:
            if row[0] != '1-1':
                continue
            if float(row[4]) / float(row[5]) < thresh or \
                    float(row[5]) / float(row[4]) < thresh:
                continue
            if s0 in third:
                lt0 = row[2].split('::')[2]
            else:
                lt0 = row[2].split('::')[1]
            if s1 in third:
                lt1 = row[3].split('::')[2]
            else:
                lt1 = row[3].split('::')[1]
            strains.add(s0); strains.add(s1)
            g.add_node(s0 + ':' + lt0)
            g.add_node(s1 + ':' + lt1)
            g.add_edge(s0 + ':' + lt0, s1 + ':' + lt1)

components = networkx.connected_components(g)
sys.stdout.write('group\t' + '\t'.join(s for s in strains) + '\n')
for i, c in enumerate(components):
    genes = collections.defaultdict(list)
    for gene in list(c):
        s = gene.split(':')[0]
        g = gene.split(':')[1]
        genes[s].append((s, g))
    rows = itertools.product(*genes.values())
    for row in rows:
        r = {r[0]:r[1] for r in row}
        sys.stdout.write(str(i))
        for s in strains:
            sys.stdout.write('\t')
            if s in r:
                sys.stdout.write(r[s])
        sys.stdout.write('\n')
