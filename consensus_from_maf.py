#!/usr/bin/env python3
"""
    Get a consensus sequence from an MAF file.

    consensus_from_maf.py [<input file>]
"""

import collections
import sys

if len(sys.argv) > 1:
    h = open(sys.argv[1], 'rt')
else:
    h = sys.stdin

with h as ih:
    annots = {}
    counters = []
    for line in ih:
        if line.startswith('a'):
            if len(counters) > 0:
                for c in counters:
                    sys.stdout.write(c.most_common()[0][0])
                sys.stdout.write('\n')
            annots = {x.split('=')[0]:x.split('=')[1] for x in line.strip().split()[1:]}
            counters = []
            sys.stdout.write('>' + annots['label'] + '\n')
        elif line.startswith('s'):
            fields = line.strip().split()
            s = fields[6]
            if len(counters) == 0:
                counters = [collections.Counter() for x in s]
            for i, b in enumerate(s):
                if b != '-':
                    counters[i].update(b)
        else:
            continue
    for c in counters:
        sys.stdout.write(c.most_common()[0][0])
    sys.stdout.write('\n')
