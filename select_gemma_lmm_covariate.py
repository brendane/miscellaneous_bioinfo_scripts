#!/usr/bin/env python2.7
"""
    Create a covariates file for the LMM mode of GEMMA using the top
    variant from a previous run. Missing genotypes are substituted by
    the mean genotype value.

    select_gemma_lmm_covariate.py --output <output file>
        --previous <previous covariates> <plink file>
        <output.assoc file>
"""

import argparse
import csv

from pyplink import PyPlink

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--previous')
parser.add_argument('plink')
parser.add_argument('assoc')
args = parser.parse_args()

top_variant = None
top_p = 1.0
with open(args.assoc, 'rb') as ihandle:
    rdr = csv.DictReader(ihandle, delimiter='\t')
    for row in rdr:
        p = float(row['p_lrt'])
        if p < top_p:
            top_p = p
            top_variant = row['rs']
if top_variant is None:
    raise Exception('No top variant chosen')
print top_variant

p = PyPlink(args.plink)
g = None
for marker, genotypes in p:
    if marker == top_variant:
        s = sum(x for x in genotypes.tolist() if x != -1)
        c = sum(1. for x in genotypes.tolist() if x != -1)
        meang = s/c
        g = [str(x) if x != -1 else str(meang) for x in genotypes]
        break
if g is None:
    raise Exception('Could not find top variant in plink file')

if args.previous is None or args.previous == 'NONE':
    previous = [['1'] for x in g]
else:
    previous = []
    with open(args.previous, 'rb') as ihandle:
        for line in ihandle:
            previous.append(line.strip().split())

with open(args.output, 'wb') as ohandle:
    for p, gg in zip(previous, g):
        ohandle.write('\t'.join(p) + '\t' + gg + '\n')
