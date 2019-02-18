#!/usr/bin/env python3
"""
    Label replicons in a fasta file by the Assigned-Molecule-Location
    column.

    label_replicon_by_assm_report.py <fasta> <report>
"""

import sys

from Bio import SeqIO

replicons = {}
with open(sys.argv[2], 'r') as ih:
    for line in ih:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        n = fields[6]
        r = fields[3]
        replicons[n] = r

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    r = replicons[rec.id]
    sys.stdout.write('>' + r + ':' + rec.id + '\n')
    sys.stdout.write(str(rec.seq) + '\n')
