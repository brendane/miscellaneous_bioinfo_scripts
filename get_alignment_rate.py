#!/usr/bin/env python2.7
import sys
for f in sys.argv[1:]:
    e = 0
    m = 0
    u = 0
    with open(f, 'r') as ih:
        for line in ih:
            if 'reads aligned to ensifer' in line:
                e = int(line.strip().split()[0])
            if 'reads aligned to medicago' in line:
                m = int(line.strip().split()[0])
            if 'reads not aligned' in line:
                u = int(line.strip().split()[0])
    print float(e + m) / (e + m + u) * 100
