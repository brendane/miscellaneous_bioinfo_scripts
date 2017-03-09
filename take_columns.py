#!/usr/bin/env python
import csv
import sys

if len(sys.argv) > 2:
    f = open(sys.argv[1])
else:
    f = sys.stdin

cols = set(sys.argv[-1].split(','))

with f:
    fieldnames = f.readline().strip().split()
    out = []
    for fn in fieldnames:
        if fn in cols:
            out.append(fn)
    print '\t'.join(out)
    for line in f:
        fields = line.strip().split()
        out = []
        for fn, ff in zip (fieldnames, fields):
            if fn in cols:
                out.append(ff)
        print '\t'.join(out)
