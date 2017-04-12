#!/usr/bin/env python2

import sys
import csv

ld_groups = dict()
with open(sys.argv[1], 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        ld_groups['group-' + row[4]] = row[3]

for line in sys.stdin:
    print ld_groups[line.strip()]
