#!/usr/bin/env python2.7
"""
    Filter a GFF file based on tags.

    gff_filter_tags.py [--only-present] <gff file> <tag1=value1,value2 ...> <tag2=...>
"""

import argparse
import csv
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--only-present', default=False, action='store_true')
parser.add_argument('gff')
parser.add_argument('filters', nargs='+')
args = parser.parse_args()

op = args.only_present
filters = {}
if args.filters[0] != "":
    filters = {x.split('=')[0]:[y for y in x.split('=')[1].split(',')]
               for x in args.filters}


with open(args.gff, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'):
            sys.stdout.write('\t'.join(row) + '\n')
            continue
        tags = {x.split('=')[0]:x.split('=')[1] for x in row[8].split(';')}
        ok = True
        for f in filters:
            if f in tags:
                if tags[f] not in filters[f]:
                    ok = False
                    break
            elif op:
                ok = False
                break
        if ok:
            sys.stdout.write('\t'.join(row) + '\n')
