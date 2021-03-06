#!/usr/bin/env python3

import csv
import sys

raise Exception('BLAT parsing script has bugs. Do not use.')

ref_pos_done = set()
qry_pos_done = set()
multi = set()
matches = set()
with open(sys.argv[1], 'rt') as ih:
    for i in range(5):
        ih.readline()
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        strand = row[8]  ## start is always < end
        q = row[9]
        r = row[13]
        qs = int(row[11])
        qe = int(row[12])
        rs = int(row[15])
        re = int(row[16])
        block_lens = [int(x) for x in row[18].split(',')[:-1]]
        block_qs = [int(x) for x in row[19].split(',')[:-1]]
        block_rs = [int(x) for x in row[20].split(',')[:-1]]
        for bl, bq, br in zip(block_lens, block_qs, block_rs):
            i = br + rs
            j = bq + qs
            if strand == '-':
                j = bq + bl - 1
            for k in range(bl):
                if (r, i) in ref_pos_done or (q, j) in qry_pos_done:
                    multi.add((r, i)); multi.add((r, j))
                    i += 1
                    if strand == '+':
                        j += 1
                    else:
                        j -= 1
                    continue
                matches.add((r, i, q, j, strand))
                ref_pos_done.add((r, i))
                qry_pos_done.add((q, j))
                i += 1
                if strand == '+':
                    j += 1
                else:
                    j -= 1

for r, i, q, j, s in matches:
    if (r, i) in multi or (q, j) in multi:
        continue
    sys.stdout.write(q + '\t' + str(j+1) + '\t' + r + '\t' + str(i+1) + '\t' + s + '\n')
