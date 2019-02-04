#!/usr/bin/env python3
"""
    Filter a fasta file so that all the entries are larger than a given
    length.

    filter_fasta_length.py <fasta input file> <threshold length>

    Writes to stdout.
"""

import sys

from Bio import SeqIO

threshold = int(sys.argv[2])
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    if len(rec) < threshold:
        continue
    SeqIO.write(rec, sys.stdout, 'fasta')
