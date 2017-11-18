#!/usr/bin/env python2.7
"""
    Remove all copies of sequences that have duplicated names in a
    fasta file.
"""

import sys

from Bio import SeqIO

names = set()
remove = set()

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    if rec.id in names:
        remove.add(rec.id)
    names.add(rec.id)

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    if rec.id in remove:
        continue
    sys.stdout.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
