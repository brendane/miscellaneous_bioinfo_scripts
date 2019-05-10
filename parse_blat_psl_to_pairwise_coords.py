#!/usr/bin/env python3

import csv
import sys

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
        #for bl, bq, br in zip(block_lens, block_qs, block_rs):
        #    if strand == '+':
        #        sys.stdout.write(q + '\t' + str(bq) + '\t' + str(bq + bl) + '\t' +
        #                         r + '\t' + str(br) + '\t' + str(br + bl) + '\n')
        #    else:
        #        sys.stdout.write(q + '\t' + str(bq + bl) + '\t' + str(bq) + '\t' +
        #                         r + '\t' + str(br) + '\t' + str(br + bl) + '\n')
        for bl, bq, br in zip(block_lens, block_qs, block_rs):
            i = br
            j = bq
            for k in range(bl):
                if (r, i) in ref_pos_done or (q, j) in qry_pos_done:
                    multi.add((r, i)); multi.add((r, j))
                    continue
                matches.add((r, i, q, j))
                ref_pos_done.add((r, i))
                qry_pos_done.add((q, j))
                #sys.stdout.write(q + '\t' + str(j+1) + '\t' + r + '\t' + str(i) + '\n')
                i += 1
                if strand == '+':
                    j += 1
                else:
                    j -= 1

for r, i, q, j in matches:
    if (r, i) in multi or (q, j) in multi:
        continue
    sys.stdout.write(q + '\t' + str(j+1) + '\t' + r + '\t' + str(i+1) + '\n')
