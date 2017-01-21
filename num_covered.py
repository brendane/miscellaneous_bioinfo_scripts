#!/usr/bin/env python2
"""
    Number of sites covered by a certain proportion of strains.

    Expects a directory with tsv files.

    prop_covered.py [options] --output <output file> <input directory>
        <degeneracy file>

        --min-depth     Minimum depth to count a site as covered (1)
        --min-prop      Minimum proportion of strains to count a site as
                        covered (0)
        --keep          File with list of strains (one per line) to
                        consider

    The degeneracy file is assumed to contain synonymous / non-synonymous
    classifications for all the sites in coding genes.
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
parser.add_argument('--min-prop', type=float, default=0.0)
parser.add_argument('--keep')
parser.add_argument('--output')
parser.add_argument('input')
parser.add_argument('degen')
args = parser.parse_args()

min_depth = args.min_depth
min_prop = args.min_prop
out_file = args.output
indir = args.input
deg_file = args.degen

site_class = {}
with open(deg_file, 'rb') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    for row in rdr:
        key = (row['replicon'], int(row['pos']))
        site_class[key] = (float(row['syn']), float(row['nonsyn']))

keep = set()
if args.keep is not None:
    with open(args.keep, 'rb')  as handle:
        for line in handle:
            keep.add(line.strip())

with open(out_file, 'wb') as out:
    out.write('feature\treplicon\tpos\tnsites\tnsites_genic\tnsites_syn\t' +
              'nsites_nonsyn\tnsites_intergenic\n')
    for fname in os.listdir(indir):
        if not fname.endswith('.tsv') and not fname.endswith('.tsv.gz'):
            continue
        feature = re.sub('\.tsv', '', osp.splitext(fname)[0])
        infile = osp.join(indir, fname)
        if fname.endswith('.gz'):
            data = pandas.read_table(gzip.open(infile, 'rb'), sep='\t', header=0)
        else:
            data = pandas.read_table(infile, sep='\t', header=0)
        strains = list(data.columns[2:])
        for strain in strains:
            if len(keep) > 0 and strain not in keep:
                data.drop(strain, axis=1, inplace=True)
        nsites = 0
        nsites_g = 0
        nsites_i = 0
        nsites_n = 0.
        nsites_s = 0.
        for row in data.iterrows():
            pcov = sum(1. for d in row[1][2:] if d >= min_depth) / \
                    (len(row[1]) - 2)
            if pcov >= min_prop:
                nsites += 1
                key = (row[1]['replicon'], row[1]['pos'])
                try:
                    s, ns = site_class[key]
                    nsites_g += 1
                    nsites_s += s
                    nsites_n += ns
                except KeyError:
                    nsites_i += 1
        out.write(feature + '\t' + data['replicon'][0] + '\t' +
                  str(data['pos'][0]) + '\t' + str(nsites) + '\t' +
                  str(nsites_g) + '\t' + str(nsites_s) + '\t' +
                  str(nsites_n) + '\t' + str(nsites_i) + '\n')
