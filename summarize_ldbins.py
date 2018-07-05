#!/usr/bin/env python2.7
"""
    summarize_ldbins.py <input file>
"""

import sys
import csv
import numpy

stats = {}

ld_bin = 0.
ld_values = []
with open(sys.argv[1], 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        lb = float(row[2])
        if lb != ld_bin:
            ld_values.sort()
            stats[ld_bin] = {}
            stats[ld_bin]['mean'] = numpy.mean(ld_values)
            stats[ld_bin]['quantiles'] = [numpy.percentile(ld_values, q) for q in range(5, 100, 5)]
            ld_bin = lb
            ld_values = []
        else:
            ld_values.append(float(row[0]))

sys.stdout.write('bin\tmean\tp05\tp10\tp15\tp20\tp25\tp30\tp35\tp40\tp45\tp50\tp55\tp60\tp65\tp70\tp75\tp80\tp85\tp90\tp95\n')
for lb in sorted(stats):
    sts = stats[lb]
    sys.stdout.write(str(lb) + '\t' + str(sts['mean']) + '\t' + 
                     '\t'.join(str(x) for x in sts['quantiles']) + '\n')
