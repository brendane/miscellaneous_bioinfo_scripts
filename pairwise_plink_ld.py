#!/usr/bin/env python2.7
"""
    Filter plink LD file

    filter_plink_ld.py input
"""

import collections
import gzip
import sys

infile = sys.argv[1]

sum_r2 = {'Uni0':0., 'Uni1':0., 'Uni2':0., 'denovo':0., 'snp_pav':0.,
          'snp_snp':0., 'pav_pav':0, 'all':0.}
count_r2 = {'Uni0':0, 'Uni1':0, 'Uni2':0, 'denovo':0, 'snp_pav':0,
            'snp_snp':0, 'pav_pav':0, 'all':0}

open_fun = open
if infile == '-':
    handle = sys.stdin
elif infile.endswith('.gz'):
    handle = gzip.open(infile, 'rb')
else:
    handle = open(infile, 'rb')
with handle as ih:
    ih.readline()
    for line in ih:
        v0, v1, r2 = line.strip().split('\t')
        r2 = float(r2)

        tp = ['all']
        if v0.startswith('snp') and v1.startswith('snp'):
            tp.append('snp_snp')
        elif v0.startswith('snp') and v1.startswith('rdv'):
            tp.append('snp_pav')
        elif v0.startswith('rdv') and v1.startswith('snp'):
            tp.append('snp_pav')
        elif v0.startswith('rdv') and v1.startswith('snp'):
            tp.append('pav_pav')
        if 'Uni0' in v0 and 'Uni0' in v1:
            tp.append('Uni0')
        elif 'Uni1' in v0 and 'Uni1' in v1:
            tp.append('Uni1')
        elif 'Uni2' in v0 and 'Uni2' in v1:
            tp.append('Uni2')
        elif 'denovo' in v0 and 'denovo' in v1:
            tp.append('denovo')

        for t in tp:
            sum_r2[t] += r2
            count_r2[t] += 1

for tp in sum_r2:
    if count_r2[tp] > 0:
        sys.stdout.write(tp + '\t' + str(sum_r2[tp] / count_r2[tp]) + '\n')
