#!/usr/bin/env python2.7
"""
    Given two files from mnp2snp.py report which variants in the
    second file match which variants in the first file.

    compare_mnp_real_coords.py <file 1> <file 2>

    If a variant in file 1 matches more than one variant in file 2,
    none of them are reported.
"""

import collections
import csv
import sys

first_file = sys.argv[1]
second_file = sys.argv[2]

first_rs = set()
first_pos = {}
with open(first_file, 'rb') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        first_rs.add(row['rs'])
        first_pos[(row['chrom'], row['pos'], row['alt_allele'])] = row['rs']

second_pos = {}
second_rs = set()
with open(second_file, 'rb') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        second_pos[(row['chrom'], row['pos'], row['alt_allele'])] = row['rs']
        second_rs.add(row['rs'])

equivalents = {}
equivalent_count = collections.defaultdict(int)
for (c, p, a), rs in second_pos.iteritems():
    if rs in first_rs:
        equivalents[rs] = rs
        equivalent_count[rs] += 1
    elif (c, p, a) in first_pos:
        rs_1 = first_pos[(c, p, a)]
        #if rs_1 in second_rs:
        #    # This variant is not the best match for its equivalent
        #    continue
        equivalents[rs] = rs_1
        equivalent_count[rs_1] += 1

for rs_2, rs_1 in equivalents.iteritems():
    #if equivalent_count[rs_1] > 1:
    #    continue
    sys.stdout.write(rs_1 + '\t' + rs_2 + '\n')
