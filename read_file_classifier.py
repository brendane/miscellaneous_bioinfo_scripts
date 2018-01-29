#!/usr/bin/env python2.7

import collections
import gzip
import os.path as osp
import re
import sys

from Bio import SeqIO

strains = {}

dups = set()
for fname in sys.argv[1:]:
    if 'R2' in fname or 'fastq' not in fname:
        continue
    name = re.sub('_.+', '', osp.basename(fname))
    ofun = open
    if fname.endswith('.gz'):
        ofun = gzip.open
    with ofun(fname, 'rb') as h:
        for rec in SeqIO.parse(h, 'fastq'):
            if ':N:' in rec.description.split()[-1]:
                first_part = re.sub(':.+', '', rec.id)
                second_part = rec.description.split()[-1]
                key = first_part + '-' + second_part
            else:
                key = re.sub('\\..+', '', rec.description.split()[0])
            if key in strains:
                dups.add(key)
            strains[key] = name
            break

for k, s in strains.iteritems():
    print k + '\t' + s
