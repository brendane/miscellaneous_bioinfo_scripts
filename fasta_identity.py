#!/usr/bin/env python2.7

import sys

from Bio import AlignIO

aln = AlignIO.read(sys.argv[1], 'fasta')

l = len(aln[0])
i = 0
l_ng = 0
for p in xrange(l):
    b0 = aln[0, p].upper()
    b1 = aln[1, p].upper()
    if b0 == b1:
        i += 1
    if b0 != '-' and b1 != '-':
        l_ng += 1

print 'ID\t' + str(i / float(l))
print 'ID ungapped\t' + str(i / float(l_ng))
