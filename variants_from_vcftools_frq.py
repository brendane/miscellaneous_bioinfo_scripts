#!/usr/bin/env python2
"""
    Get variant positions from a vcftools frq file.

    variants_from_vcftools_frq.py <input file> <output file>
"""

#==============================================================================#

import csv
import sys

#==============================================================================#


with open(sys.argv[2], 'wb') as out:
    with open(sys.argv[1], 'rb') as inhandle:
        rdr = csv.reader(inhandle, delimiter='\t')
        rdr.next()
        for row in rdr:
            freqs = [float(x.split(':')[1]) for x in row[4:]]
            if max(freqs) == 1.0:
                continue
            out.write(row[0] + '\t' + str(int(row[1])-1) + '\t' + row[1]
                      + '\t' + row[0] + '-' + row[1] + '\n')
