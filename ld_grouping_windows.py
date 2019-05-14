#!/usr/bin/env python3
"""
    Calculate proportion of shared SNPs in LD across genomic windows.

    ld_grouping_windows.py <ld file> <window size>

    The ld file is the output from r2_groups_xtra_sort, sorted by
    group.
"""

import collections
import csv
import itertools
import sys

ld_file = sys.argv[1]
win_size = int(sys.argv[2])

snp_count = collections.defaultdict(int)
shared_count = collections.defaultdict(int)
done = set()
with open(ld_file, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for group, rows in itertools.groupby(rdr, lambda x: x['group']):
        rows = list(rows)
        for r1, r2 in itertools.combinations(rows, 2):
            rs1 = r1['rs']
            rs2 = r2['rs']
            w1 = (r1['chrom'], int(r1['pos']) // win_size)
            w2 = (r2['chrom'], int(r2['pos']) // win_size)
            if rs1 not in done:
                snp_count[w1] += 1
                done.add(rs1)
            if rs2 not in done:
                snp_count[w2] += 1
                done.add(rs2)
            if w1[1] > w2[1]:
                shared_count[(w1, w2)] += 1
            else:
                shared_count[(w2, w1)] += 1

for (w1, w2), c in shared_count.items():
    ## Proportion of possible SNP pairs that are in LD
    if w1 == w2:
        p =  c / (snp_count[w1] * (snp_count[w1]) / 2)
    else:
        p =  c / (snp_count[w1] * snp_count[w2])
    sys.stdout.writelines([w1[0], '\t', str(w1[1]), '\t',
                           w2[0], '\t', str(w2[1]), '\t',
                           str(p), '\t', str(c), '\n'])
