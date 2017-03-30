#!/usr/bin/env python2.7

import sys
infiles = sys.argv[1:]
data = []
genes = []
for i, fname in enumerate(infiles):
    sys.stderr.write(fname + '\n')
    d = []
    with open(fname, 'rb') as ihandle:
        for j, line in enumerate(ihandle):
            g, c = line.strip().split()
            if i != 0 and g != genes[j]:
                raise Exception('no match')
            if i == 0:
                genes.append(g)
            d.append(c)
    data.append(d)
out = sys.stdout
for i in xrange(len(genes)):
    out.write(genes[i])
    for d in data:
        out.write('\t' + d[i])
    out.write('\n')
