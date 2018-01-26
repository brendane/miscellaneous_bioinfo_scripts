#!/usr/bin/env python2.7
"""
    Filter a DGRP file.

    filter_dgrp.py --output <output file> --min-gt <minimum gt rate>
        <input file>
"""

import argparse
import csv

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--min-gt', type=float, default=0.)
parser.add_argument('infile')
args = parser.parse_args()

with open(args.output, 'wb') as oh:
    with open(args.infile, 'rb') as ih:
        oh.write(ih.readline())
        rdr = csv.reader(ih)
        for row in rdr:
            n = len(row) - 2
            ngt = sum(1. for b in row[2:] if b.upper() != 'N')
            if ngt / n >= args.min_gt:
                oh.write(','.join(row) + '\n')
