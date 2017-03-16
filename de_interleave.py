#!/usr/bin/env python2.7

import itertools
import sys
import re

from Bio import SeqIO

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

with open(sys.argv[1] + '_R1.fastq', 'wb') as out1:
    with open(sys.argv[1] + '_R2.fastq', 'wb') as out2:
        for rec1, rec2 in grouper(SeqIO.parse(sys.stdin, 'fastq'), 2):
            if re.sub('/1$', '', rec1.name) != re.sub('/2$', '', rec2.name):
                raise Exception(rec1.name + ' ' + rec2.name + ' do not match')
            SeqIO.write(rec1, out1, 'fastq')
            SeqIO.write(rec2, out2, 'fastq')
