#!/usr/bin/env python2.7

import csv
import itertools
import sys

with open(sys.argv[1], 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for qname, rows in itertools.groupby(rdr, lambda x: x[0]):
        top_hit = rows.next()
        print top_hit[0] + '\t' + ' '.join(top_hit[5].split()[:2])
