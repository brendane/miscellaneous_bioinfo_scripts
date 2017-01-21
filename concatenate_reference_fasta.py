#!/usr/bin/env python
"""
    Concatenate the records in a FASTA file with 100 Ns in
    between.

    concatenate_reference_fasta.py <infile> <contig order>

"""

#==============================================================================#

from sys import stdout, argv
import re

from Bio import SeqIO

#==============================================================================#

gap = 100

infile = argv[1]
idx = SeqIO.index(infile, 'fasta')

sq_order = argv[2:]
sq_lens = {}
padding = {}
for id, rec in idx.iteritems():
    sq_lens[id] = len(rec)
p = 0
for sq in sq_order:
    padding[sq] = p
    p += sq_lens[sq] + gap

stdout.write('>genome\n')
first = True
for id in sq_order:
    rec = idx[id]
    if first:
        first = False
    else:
        stdout.write('N' * gap)
    stdout.write(str(rec.seq))
stdout.write('\n')
