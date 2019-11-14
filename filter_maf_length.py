#!/usr/bin/env python3
"""
    Filter an MAF file based on LCB length.

    consensus_from_maf.py <min length> <input file>
"""

import collections
import sys

min_len = int(sys.argv[1])

if len(sys.argv) > 2:
    h = open(sys.argv[2], 'rt')
else:
    h = sys.stdin

with h as ih:
    annots = {}
    keep = True
    data = []
    for line in ih:
        if line.startswith('a'):
            if len(data) > 0 and keep:
                for l in data:
                    sys.stdout.write(l)
                sys.stdout.write('\n')
            data.append(line)
        elif line.startswith('s'):
            if len(line.strip().split()[6]) < min_len:
                keep = False
                continue
            data.append(line)
    if len(data) > 0 and keep:
        for l in data:
            sys.stdout.write(l)
