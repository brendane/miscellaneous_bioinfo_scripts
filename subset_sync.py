#!/usr/bin/env python2
"""
    Use a bed file to subset a Popoolation synchronized file.

    subset_sync.py <bed file> <output file> <sync files...>
"""

#==============================================================================#

import csv
import gzip
import sys

#==============================================================================#

print 'Reading positions file'
positions = set()
with open(sys.argv[1], 'rb') as handle:
    for line in handle:
        contig, _, pos = line.strip().split()[:3]
        positions.add((contig, pos))

print 'Subsetting'
with open(sys.argv[2], 'wb') as out:
    for fname in sys.argv[3:]:
        opfun = open
        if fname.endswith('.gz'):
            opfun = gzip.open
        with opfun(fname, 'rb') as handle:
            rdr = csv.reader(handle, delimiter='\t')
            for row in rdr:
                if (row[0], row[1]) not in positions:
                    continue
                out.write('\t'.join(row) + '\n')
