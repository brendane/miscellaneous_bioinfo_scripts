#!/usr/bin/env python2.7
"""
    Filter a fasta file by missingness. Does not count gaps as missing.

    fasta_missingness_filter.py <infile> <max missing>

    Sites with too much missing data are not written.
"""

import sys

from Bio import AlignIO

mm = float(sys.argv[2])

aln = AlignIO.read(sys.argv[1], 'fasta')
n_strains = float(len(aln))

missing = set()
for i in xrange(aln.get_alignment_length()):
    if (sum(1 if b.upper() == 'N' else 0 for b in aln[:, i]) / n_strains) > mm:
        missing.add(i)

for rec in aln:
    sys.stdout.write('>' + rec.id + '\n')
    for i in xrange(len(rec)):
        if i in missing:
            continue
        else:
            sys.stdout.write(rec[i])
    sys.stdout.write('\n')
