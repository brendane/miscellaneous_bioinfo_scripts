#!/usr/bin/env python2
"""
    Proportion covered.

    Expects a directory with tsv files.

    prop_covered.py [options] --output <output file> <input directory>

        --min-depth     Minimum depth to count a site as covered (1)
"""

import argparse
import csv
import gzip
import os
import os.path as osp
import re

import pandas

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--min-depth', type=int, default=1)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

min_depth = args.min_depth
out_file = args.output
indir = args.input

covered = {}
for fname in os.listdir(indir):
    if not fname.endswith('.tsv') and not fname.endswith('.tsv.gz'):
        continue
    feature = re.sub('\.tsv', '', osp.splitext(fname)[0])
    infile = osp.join(indir, fname)
    if fname.endswith('.gz'):
        data = pandas.read_table(gzip.open(infile, 'rb'), sep='\t', header=0)
    else:
        data = pandas.read_table(infile, sep='\t', header=0)
    covd = {}
    strains = list(data.columns[2:])
    for strain in strains:
        covd[strain] = sum(1 for d in data[strain] if d >= min_depth) / float(len(data[strain]))
    covered[(feature, data['replicon'][0], data['pos'][0])] = covd

with open(out_file, 'wb') as out:
    out.write('feature\treplicon\tpos\t' + '\t'.join(strains) + '\n')
    for feature in covered:
        out.write(feature[0] + '\t' + feature[1] + '\t' + str(feature[2]))
        for strain in strains:
            c = covered[feature][strain]
            out.write('\t' + str(c))
        out.write('\n')
