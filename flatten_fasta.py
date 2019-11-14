#!/usr/bin/env python3

import sys

from Bio import SeqIO

for rec in SeqIO.parse(sys.stdin, 'fasta'):
    print('>', rec.description)
    print(str(rec.seq))
