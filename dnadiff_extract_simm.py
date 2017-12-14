#!/usr/bin/env python2.7
"""
    Extract similarity information from a dnadiff report:
    - Percent similarity in 1-1 alignment regions
    - Percent similarity in M-M alignment regions
    - Proportion of shared genome =
        (aligned bases 1 + aligned bases 2) /
        (genome size 1 + genome size 2)
"""

import re
import sys

bases = False
oo = False
mm = False
n_total_bases = 0
n_total_aligned_bases = 0
prop_shared = 0.
oo_simm = 0.
mm_simm = 0.
with open(sys.argv[1], 'rb') as ih:
    for line in ih:
        line = line.strip()
        if bases:
            if line.startswith('TotalBases'):
                n0, n1 = map(int, line.split()[1:])
                n_total_bases = n0 + n1
            elif line.startswith('AlignedBases'):
                a0, a1 = map(lambda x: float(re.sub('\(.+', '', x)),
                             line.split()[1:])
                prop_shared = (a0 + a1) / n_total_bases
        elif oo:
            if line.startswith('AvgIdentity'):
                oo_simm = float(line.split()[1])
        elif mm:
            if line.startswith('AvgIdentity'):
                mm_simm = float(line.split()[1])
        if line == '[Bases]':
            bases = True
            oo = False
            mm = False
        elif line.startswith('1-to-1'):
            bases = False
            oo = True
            mm = False
        elif line.startswith('M-to-M'):
            bases = False
            oo = False
            mm = True

sys.stdout.write(str(oo_simm) + '\t' + str(mm_simm) + '\t' +
                 str(round(prop_shared, 4)*100))
