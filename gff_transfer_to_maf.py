#!/usr/bin/env python3
"""
    Lift over annotation coordinates to a whole genome alignment.

    gff_transfer_to_maf.py [--sep <separator (.)>] <gff file> <input file> <strain>
"""

import argparse
import collections
import csv

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--sep', default='.')
parser.add_argument('gff')
parser.add_argument('maf')
parser.add_argument('strain')
args = parser.parse_args()

coords = {} # Not sure the best way to store coordinates;
            # maybe key by original contig, then list sorted by start
            # in original.
with open(args.maf, 'rt') as ih:
    for line in ih:
        if line.startswith('a'):
            annots = {x.split('=')[0]:x.split('=')[1] for x in line.strip().split()[1:]}
        elif line.startswith('s'):
            fields = line.strip().split()
            s = fields[6]
        else:
            continue
