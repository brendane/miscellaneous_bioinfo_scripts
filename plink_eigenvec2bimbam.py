#!/usr/bin/env python2.7
"""
    Convert a plink eigenvec output file (in tab-delimited format) into
    a bimbam file for GEMMA.

    plink_eigenvec2bimbam.py <input eigenvec> <input phenotype> 
        <output prefix> <phenotype> <n axes>

    The phenotype file should have a column labeled "strain" and should
    be tab-delimited.
"""

import csv
import sys

n = int(sys.argv[5])
pheno = sys.argv[4]

data = {}
strains = []
with open(sys.argv[1], 'rb') as ihnd:
    for line in ihnd:
        fields = line.strip().split('\t')
        strains.append(fields[1])
        data[fields[1]] = ([float(f) for f in fields[2:]])

pheno_data = {}
with open(sys.argv[2], 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    for row in rdr:
        if row[pheno] == '' or 'NA' in row[pheno] or \
           'nan' in row[pheno] or 'Na' in row[pheno]:
            continue
        pheno_data[row['strain']] = row[pheno]

strains.sort()
with open(sys.argv[3] + '.geno', 'wb') as ohnd:
    for i in range(n):
        ohnd.write('pc' + str(i) + ', A, T')
        vals = [data[strain][i] for strain in strains]
        range_vals = max(vals) - min(vals)
        norm_vals = [(v - min(vals)) / range_vals * 2 for v in vals]
        for j, strain in enumerate(strains):
            if strain not in pheno_data:
                continue
            ohnd.write(', ' + str(norm_vals[j]))
        ohnd.write('\n')

with open(sys.argv[3] + '.pheno', 'wb') as ohnd:
    for strain in strains:
        if strain not in pheno_data:
            continue
        ohnd.write(pheno_data[strain] + '\n')

with open(sys.argv[3] + '.strains', 'wb') as ohnd:
    for strain in strains:
        if strain not in pheno_data:
            continue
        ohnd.write(strain + '\n')
