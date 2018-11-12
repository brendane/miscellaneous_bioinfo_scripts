#!/usr/bin/env python2.7
"""
    Filter a spades assembly by length. Renames contigs to ctg1,
    ctg2, etc.

    spades_length_filter.py <input file> <min length>
"""

import sys

from Bio import SeqIO

m = int(sys.argv[2])

for i, rec in enumerate(SeqIO.parse(sys.argv[1], 'fasta')):
    if len(rec) < m:
        continue
    sys.stdout.write('>ctg' + str(i+1) + '\n')
    sys.stdout.write(str(rec.seq) + '\n')
