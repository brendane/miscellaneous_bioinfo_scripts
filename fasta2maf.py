#!/usr/bin/env python3

import os.path as osp
import sys

from Bio import SeqIO

for i, rec in enumerate(SeqIO.parse(sys.argv[1], 'fasta')):
    sys.stdout.write('a label=' + str(i) + '\n')
    sys.stdout.write('s\t' + sys.argv[2] + '.' + rec.id + '\t' +
                     '0\t' + str(len(rec)) + '\t+\t' +
                     str(len(rec)) + str(rec.seq) + '\n')
