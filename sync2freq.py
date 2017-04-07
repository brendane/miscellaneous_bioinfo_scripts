#!/usr/bin/env python
"""
    Convert a sync file and associated list of pools into a table of
    reference allele frequencies. Excludes ambiguous reads.

    sync2freq.py <sync file> <pool file>

    Output goes to stdout.
"""

import csv
import sys

sync_file = sys.argv[1]
pool_file = sys.argv[2]

pools = []
with open(pool_file, 'rb') as ih:
    for line in ih:
        pools.append(line.strip())

oh = sys.stdout
oh.write('\t'.join(pools) + '\n')
with open(sync_file, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        first = True
        ref = row[2]
        for pool, counts in zip(pools, row[3:]):
            count_dict = {b:c for b, c in zip(['A','T','C','G','N','D'],
                                              map(int, counts.strip().split(':')))}
            ref_count = count_dict[ref]
            total_count = sum(count_dict.itervalues()) - count_dict[ref] - \
                    count_dict['N'] + ref_count
            if total_count == 0:
                ref_freq = float('nan')
            else:
                ref_freq = float(ref_count) / total_count
            if not first:
                oh.write('\t')
            first = False
            oh.write(str(ref_freq))
        oh.write('\n')
