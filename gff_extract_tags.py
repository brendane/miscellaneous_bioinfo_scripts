#!/usr/bin/env python

import csv
import sys

fields = sys.argv[2].split(',')

with open(sys.argv[1], 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'):
            sys.stdout.write('\t'.join(row) + '\n')
            continue
        tags = {x.split('=')[0]:x.split('=')[1] for x in row[8].split(';')}
        sys.stdout.write('\t'.join(row[:8]) + '\t')
        t = []
        for f in fields:
            if f in tags:
                t.append(f + '=' + tags[f] + ';')
            if ':' + f.lower() in tags:
                t.append(f + '=' + tags[':' + f.lower()] + ';')
        sys.stdout.write(''.join(t) + '\n')
