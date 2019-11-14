#!/usr/bin/env python3
"""
    Filter an MAF file based on LCB length.

    consensus_from_maf.py <min length> [<input file>]
"""

import argparse
import collections
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--relabel-by-number', default=False, action='store_true')
parser.add_argument('argv', nargs='+')
args = parser.parse_args()

min_len = int(args.argv[0])

if len(args.argv) > 1:
    h = open(args.argv[1], 'rt')
else:
    h = sys.stdin

with h as ih:
    annots = {}
    keep = True
    data = []
    label_num = 0
    for line in ih:
        if line.startswith('a'):
            if len(data) > 0 and keep:
                if args.relabel_by_number:
                    sys.stdout.write('a label=%i\n' % label_num)
                else:
                    sys.stdout.write(data[0])
                for l in data[1:]:
                    sys.stdout.write(l)
                sys.stdout.write('\n')
                label_num += 1
            data = []
            keep = True
            data.append(line)
        elif line.startswith('s'):
            if len(line.strip().split()[6]) < min_len:
                keep = False
            else:
                data.append(line)
    if len(data) > 0 and keep:
        for l in data:
            sys.stdout.write(l)
