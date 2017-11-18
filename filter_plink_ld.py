#!/usr/bin/env python2.7
"""
    Filter plink LD file

    filter_plink_ld.py min_cluster_size min_r2_with_others
        min_other_cluster input file one_var_file
"""

import collections
import gzip
import sys

min_cl_size = int(sys.argv[1])
min_r2 = float(sys.argv[2])
min_clusts = int(sys.argv[3])
infile = sys.argv[4]
ov_file = sys.argv[5]

sys.stderr.write('\treading one variant file\n')
seeds = {}
with open(ov_file, 'rb') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        seed = fields[3]
        n_vars = len(fields[5].split(','))
        seeds[seed] = n_vars

n_connections = collections.defaultdict(int)
ld_data = []

sum_r2 = {'Uni0':0., 'Uni1':0., 'Uni2':0., 'denovo':0., 'snp_pav':0.,
          'snp_snp':0., 'pav_pav':0, 'all':0.}
count_r2 = {'Uni0':0, 'Uni1':0, 'Uni2':0, 'denovo':0, 'snp_pav':0,
            'snp_snp':0., 'pav_pav':0, 'all':0.}

k = 0
sys.stderr.write('\treading LD file\n')
open_fun = open
if infile.endswith('.gz'):
    open_fun = gzip.open
with open_fun(infile, 'rb') as ih:
    ih.readline()
    for line in ih:
        v0, v1, r2 = line.strip().split('\t')
        r2 = float(r2)
        k += 1
        if k % 1000000 == 0:
            sys.stderr.write(str(len(n_connections)) + '\t')
            sys.stderr.write(str(k) + '\n')

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

        if v0 in seeds and v1 in seeds:
            #ld_data.append((v0, v1, r2))
            if r2 >= min_r2:
                n_connections[v0] += 1
                n_connections[v1] += 1
 
sys.stderr.write('\tfinding seeds to keep\n')
keep_seeds = set()
for seed, group_size in seeds.iteritems():
    if group_size < min_cl_size and n_connections[seed] < min_clusts:
        continue
    keep_seeds.add(seed)

sys.stderr.write('\twriting output\n')
with open_fun(infile, 'rb') as ih:
    ih.readline()
    for line in ih:
        v0, v1, r2 = line.strip().split('\t')
        if v0 in keep_seeds and v1 in keep_seeds:
            sys.stdout.write(v0 + '\t' + v1 + '\t' + r2 + '\n')

for tp in sum_r2:
    if count_r2[tp] > 0:
        sys.stderr.write(tp + '\t' + str(sum_r2[tp] / count_r2[tp]) + '\n')
sys.stderr.write('n seeds kept\t' + str(len(keep_seeds)) + '\n')
