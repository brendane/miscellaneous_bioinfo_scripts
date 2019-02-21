#!/usr/bin/env python3
"""
    Extract top hit (by bitscore) from human-readable NCBI BLAST output.
    Works with BLAST v2.8.1+

    blast_fmt1_top_hits.py <input file>
"""

import sys

query = None
in_alns = False
started_alns = False
hits = []
sys.stdout.write('query\thit\tscore\tevalue\n')
with open(sys.argv[1], 'r') as ih:
    for line in ih:
        if in_alns and not started_alns and line.strip() == '':
            started_alns = True
        elif in_alns and started_alns and line.strip() == '':
            in_alns = False
            started_alns = False
        elif in_alns:
            hit = line.strip().split(' | ')[0]
            score = float(line.strip().split()[-2])
            evalue = float(line.strip().split()[-1])
            hits.append((score, evalue, hit))
        if line.startswith('Query= '):
            if query is not None:
                hits.sort(key=lambda x: x[0], reverse=True)
                if len(hits) > 0:
                    sys.stdout.write(query + '\t' + hits[0][2] + '\t' +
                                     str(hits[0][0]) + '\t' +
                                     str(hits[0][1]) + '\n')
            query = ' '.join(line.strip().split()[1:])
            hits = []
        if line.startswith('Sequences producing significant alignments:'):
            in_alns = True
    hits.sort(key=lambda x: x[0])
    sys.stdout.write(query + '\t' + hits[0][2] + '\t' +
                     str(hits[0][0]) + '\t' +
                     str(hits[0][1]) + '\n')
