#!/usr/bin/env python2
"""
    Convert GFF format to bed format.

    gff2bed.py [options] <input file(s)>
        
        --include   Comma-delimited list of feature types to include
        --exclude   Comma-delimited list of feature types to exclude
        --id        Extract just the id field
"""

import argparse
import csv
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--include')
parser.add_argument('--exclude')
parser.add_argument('--id', default=False, action='store_true')
parser.add_argument('input', nargs='+')
args = parser.parse_args()

if args.include is None:
    included = None
else:
    included = set(args.include.strip().split(','))

if args.exclude is None:
    excluded = None
else:
    excluded = set(args.exclude.strip().split(','))

extract_id = args.id

for infile in args.input:
    with open(infile, 'rb') as handle:
        rdr = csv.reader(handle, delimiter='\t')
        for row in rdr:
            if row[0].startswith('#'):
                continue
            contig = row[0]
            feature_type = row[2]
            start = int(row[3]) - 1
            end = int(row[4])
            name = row[8]
            if extract_id:
                name = [x for x in name.split(';') if x.startswith('ID=')][0]
                name = re.sub('ID=', '', name)
            if included is not None and feature_type not in included:
                continue
            if excluded is not None and feature_type in excluded:
                continue
            sys.stdout.write(contig + '\t' + str(start) + '\t' + str(end) +
                             '\t' + name + '\n')
