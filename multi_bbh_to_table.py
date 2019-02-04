#!/usr/bin/env python3
"""
    1. Filter out 1-M and M-M orthologs
    2. Filter out orthologs that differ by more than 5% in length
    3. Extract locus tags

    Note that the final set of genes may include orthologs that are
    not 1-1 in two ways:
    - Strain A gene 1 and strain B gene 1 are BBHs, strain A gene 2
      and strain C gene 2 are BBHs, and strain B gene 2 and strain C
      gene 2 are BBHs. Then you get a group with A1, B1, A2, B2, C2.
    - A1, A2, B1 are in a 1-M group when analyzed pairwise, but
      B1 - C1, and A1 - C1 are 1-1. Then A1, B1, and C1 will all
      be in a group in the output, and A2 will be excluded.
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
lengths = {}
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
            lengths[(s0, lt0)] = row[4]
            lengths[(s1, lt1)] = row[5]
            strains.add(s0); strains.add(s1)
            g.add_node(s0 + ':' + lt0)
            g.add_node(s1 + ':' + lt1)
            g.add_edge(s0 + ':' + lt0, s1 + ':' + lt1)

components = networkx.connected_components(g)
sys.stdout.write('group\t' + '\t'.join(s for s in strains) + '\t' +
    '\t'.join(s + '_length' for s in strains) + '\n')
for i, c in enumerate(components):
    genes = collections.defaultdict(list)
    for gene in list(c):
        s = gene.split(':')[0]
        g = gene.split(':')[1]
        genes[s].append((s, g))
        if len(genes[s]) > 1:
            sys.stderr.write('WARNING: group %i is not 1-1\n' % i)
    rows = itertools.product(*genes.values())
    for row in rows:
        r = {r[0]:r[1] for r in row}
        l = {r[0]:lengths[r] for r in row}
        sys.stdout.write(str(i))
        for s in strains:
            sys.stdout.write('\t')
            if s in r:
                sys.stdout.write(r[s])
        for s in strains:
            sys.stdout.write('\t')
            if s in r:
                sys.stdout.write(l[s])
        sys.stdout.write('\n')
