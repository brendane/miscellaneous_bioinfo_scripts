#!/usr/bin/env python3
"""
    Produce a table of replicon assignments from show-coords output
    with the -q -T and -l options. Assumes that the reference contigs
    are named like so: <replicon>:<whatever>.

    It also assumes that the alignments were filtered with delta-filter
    -q -o 0 (these are both important).

    mummer_replicon_assignments.py <input file>
"""

#==============================================================================#

import collections
import csv
import itertools
import re
import sys

#==============================================================================#

replicons = set()
lengths = {}
cov = {}
with open(sys.argv[1], 'r') as handle:
    for i in range(4): handle.readline()
    rdr = csv.reader(handle, delimiter='\t')
    for qry, rows in itertools.groupby(rdr, lambda x: x[10]):
        qry_len = None
        qry_cov = collections.defaultdict(int)
        for row in rows:
            ref = re.sub(':.+', '', row[9])
            replicons.add(ref)
            qry_cov[ref] += int(row[5])
            ql = float(row[8])
            if qry_len is None:
                qry_len = ql
            elif ql != qry_len:
                raise Exception('lengths do not match for %s' % qry)
        lengths[qry] = qry_len
        cov[qry] = qry_cov
replicons = sorted(replicons)


sys.stdout.write('contig\t' + '\t'.join(replicons) + '\tother_unassigned\n')
for qry in cov:
    qry_cov = cov[qry]
    qry_len = lengths[qry]
    sys.stdout.write(qry)
    c = 0.
    for r in replicons:
        cv = qry_cov[r] / qry_len
        c += cv
        sys.stdout.write('\t' + str(cv))
    ou = 1 - c
    sys.stdout.write('\t' + str(ou) + '\n')
