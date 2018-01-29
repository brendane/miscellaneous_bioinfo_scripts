#!/usr/bin/env python2.7

import collections
import re
import sys

from Bio import SeqIO

counts = collections.defaultdict(int)

if len(sys.argv) > 2:
    read_keys = {}
    with open(sys.argv[2], 'rb') as ih:
        for line in ih:
            k, s = line.strip().split('\t')
            if k not in read_keys:
                read_keys[k] = [s]
            else:
                read_keys[k].append(s)

total = 0.
for rec in SeqIO.parse(sys.argv[1], 'fastq'):
    if ':N:' in rec.description.split()[-1]:
        first_part = re.sub(':.+', '', rec.id)
        second_part = rec.description.split()[-1]
        key = first_part + '-' + second_part
        if len(sys.argv) > 2:
            ss = read_keys[key]
            for s in ss:
                counts[s] += 1
        else:
            counts[key] += 1
    else:
        key = re.sub('\\..+', '', rec.description.split()[0])
        if len(sys.argv) > 2:
            ss = read_keys[re.sub('\\..+', '', rec.description.split()[0])]
            for s in ss:
                counts[s] += 1
        else:
            counts[key] += 1
    total += 1.

for k, v in counts.iteritems():
    print k + '\t' + str(v) + '\t' + str(v / total)
