#!/usr/bin/env python2.7
"""
    Convert strain identity file from dnadiff_extract_simm.py to
    phylip-style distance matrix. Subtracts all values from 100 to get
    a true distance matrix.
"""

import csv
import sys

keep = set()
if len(sys.argv) > 2:
    with open(sys.argv[2], 'rb') as ih:
        for line in ih:
            keep.add(line.strip())

dists = {}
with open(sys.argv[1], 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        s1, s2 = row[:2]
        if len(keep) > 0 and (s1 not in keep or s2 not in keep):
            continue
        dist = (100 - float(row[2])) / 100.
        if s1 not in dists:
            dists[s1] = {}
        dists[s1][s2] = dist
        if s2 not in dists:
            dists[s2] = {}
        dists[s2][s1] = dist

strains = sorted(dists)
sys.stdout.write(str(len(strains)) + '\n')
for s1 in strains:
    sys.stdout.write(s1)
    for s2 in strains:
        if s1 == s2:
            d = '0'
        else:
            d = str(dists[s1][s2])
        sys.stdout.write('\t' + d)
    sys.stdout.write('\n')
