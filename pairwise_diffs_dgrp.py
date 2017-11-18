#!/usr/bin/env python2.7
"""
    Count pairwise differences in a haploid DGRP format file.

    pairwise_diffs_dgrp.py <input file>
"""

import collections
import csv
import sys

pairwise_diffs = collections.defaultdict(int)
with open(sys.argv[1], 'rb') as ih:
    rdr = csv.DictReader(ih, delimiter=',')
    strains = rdr.fieldnames[2:]
    for row in rdr:
        for i, s0 in enumerate(strains[:-1]):
            for s1 in strains[(i+1):]:
                if row[s0].upper() == 'N' or row[s1].upper() == 'N':
                    continue
                if row[s0].upper() != row[s1].upper():
                    pairwise_diffs[(s0, s1)] += 1

for i, s0 in enumerate(strains[:-1]):
    for s1 in strains[(i+1):]:
        print '%s\t%s\t%i' % (s0, s1, pairwise_diffs[(s0, s1)])
