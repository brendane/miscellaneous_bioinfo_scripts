#!/usr/bin/env python2.7
"""
    Merge GEMMA BSLMM runs.
"""

import csv
import os
import os.path as osp
import re
import sys

iodir = sys.argv[1]

hyp_files = []
param_files = []
for fname in os.listdir(iodir):
    if fname.endswith('param.txt'):
        param_files.append(fname)
    if fname.endswith('hyp.fixed.txt'):
        hyp_files.append(fname)

pdata = {}
rs = []
for i, fname in enumerate(param_files):
    run = re.sub('\\..+', '', re.sub('output_', '', fname))
    with open(osp.join(iodir, fname), 'rb') as ihandle:
        rdr = csv.DictReader(ihandle, delimiter='\t')
        for row in rdr:
            if row['rs'] not in pdata:
                pdata[row['rs']] = []
            pdata[row['rs']].append(row)
            if i == 0:
                rs.append(row['rs'])

with open(osp.join(iodir, 'output.param.merged.tsv'), 'wb') as out:
    out.write('chr\trs\tps\tn_miss\talpha\tbeta\tgamma\talphas\tbetas\tgammas\n')
    for r in rs:
        r_data = pdata[r]
        alphas = [float(d['alpha']) for d in r_data]
        betas = [float(d['beta']) for d in r_data]
        gammas = [float(d['gamma']) for d in r_data]
        mean_alpha = sum(alphas) / len(alphas)
        mean_beta = sum(betas) / len(betas)
        mean_gamma = sum(gammas) / len(gammas)
        out.write(r_data[0]['chr'] + '\t' +
                  r + '\t' +
                  r_data[0]['ps'] + '\t' +
                  r_data[0]['n_miss'] + '\t' +
                  str(mean_alpha) + '\t' +
                  str(mean_beta) + '\t' +
                  str(mean_gamma) + '\t' +
                  ','.join(str(f) for f in alphas) + '\t' +
                  ','.join(str(f) for f in betas) + '\t' +
                  ','.join(str(f) for f in gammas) +
                  '\n')

with open(osp.join(iodir, 'output.hyp.merged.tsv'), 'wb') as out:
    for i, fname in enumerate(hyp_files):
        with open(osp.join(iodir, fname), 'rb') as ihandle:
            line = ihandle.readline()
            if i == 0:
                out.write('run\t' + line.strip() + '\n')
            for line in ihandle:
                out.write(str(i) + '\t' + line.strip() + '\n')
