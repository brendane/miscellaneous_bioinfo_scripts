#!/usr/bin/env python2.7
"""
    Convert a gff3 file into a tsv file of annotation information.

    gff_annotation_table.py [--strain <strain name>] <infile>
        <column list>
"""

import argparse
import csv
import re
import sys

def gff_info(info_col):
    return {x.split('=')[0]:x.split('=')[1]
            for x in info_col.strip().split(';')}

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--strain')
parser.add_argument('infile')
parser.add_argument('columns')
args = parser.parse_args()

strain = args.strain
cols = args.columns.split(',')

if strain is not None:
    sys.stdout.write('strain\t')
sys.stdout.write('\t'.join(cols) + '\n')
with open(args.infile, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'):
            continue
        data = {}
        data['contig'] = row[0]
        data['start'] = row[3]
        data['end'] = row[4]
        data['strand'] = row[6]
        data.update(gff_info(row[8]))
        data['annotation'] = ''
        if 'gene' in data:
            data['annotation'] += data['gene']
        if 'product' in data:
            data['annotation'] += ': ' + data['product']
        if strain is not None:
            sys.stdout.write(strain + '\t')
        output = []
        for c in cols:
            if c in data:
                output.append(data[c])
            else:
                output.append('')
        sys.stdout.write('\t'.join(output) + '\n')
