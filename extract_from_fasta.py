#!/usr/bin/env python3
"""
    Extract fasta records given a file with the record names.

    extract_from_fasta.py <text file with IDs> <fasta file(s)>
"""

import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('ids')
parser.add_argument('fastas', nargs='+')
args = parser.parse_args()

ids = set()
with open(args.ids, 'r') as ih:
    for line in ih:
        ids.add(line.strip())

for fname in args.fastas:
    for rec in SeqIO.parse(fname, 'fasta'):
        if rec.description in ids:
            print('>' + rec.description)
            print(str(rec.seq))
