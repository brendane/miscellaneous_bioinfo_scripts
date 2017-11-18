#!/usr/bin/env python2.7
"""
    Given a three column tsv file, extract BBHs. Expects taxon name to
    be the first element separated by "::".

    extract_bbh.py <input file>
"""

import csv
import itertools
import sys

import networkx as nx

# Read in the data
scores = {}
seqs = set()
with open(sys.argv[1], 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        s0, s1 = row[:2]
        score = float(row[2])
        seqs.add(s0)
        seqs.add(s1)
        if s0 not in scores:
            scores[s0] = {}
        scores[s0][s1] = score

# Get the best hit for each sequence
best_hit_graph = nx.Graph()
best_hits = {}
for s in seqs:
    try:
        hits = scores[s]
    except KeyError:
        sys.stderr.write('Skipping %s\n' %s)
    best_score = max(hits.values())
    best_hits[s] = []
    for s1, scr in hits.iteritems():
        if scr == best_score:
            best_hit_graph.add_edge(s, s1)
            best_hits[s].append((s1, best_score))

# Sort into 1-1, 1-many, and many-many groups
#
# 1-1: Each sequence is the other's best hit and is not the best hit for
#      any other sequence
# 1-M: One sequence is the best hit for several sequences in the other
#      sample; could be a duplication or a deletion
# M-M: Same as above, but in both directions or a more complicated
#      scenario.
components = nx.connected_components(best_hit_graph)
o_o = []
o_m = []
m_m = []
o_o_seqs = set()
o_m_seqs = set()
m_m_seqs = set()
taxa = set()
for c in components:
    if len(taxa) == 0:
        for s in c:
            t = s.split('::')[0]
            taxa.add(t)
    if len(c) == 2:
        o_o.append(c)
        for s in c:
            o_o_seqs.add(s)
    else:
        taxon_count = {}
        for s in c:
            t = s.split('::')[0]
            if t not in taxon_count:
                taxon_count[t] = 0
            taxon_count[t] += 1
        if min(taxon_count.values()) == 1:
            o_m.append(c)
            for s in c:
                o_m_seqs.add(s)
        else:
            m_m.append(c)
            for s in c:
                m_m_seqs.add(s)

taxa = sorted(taxa)
if len(taxa) != 2:
    raise Exception('Wrong number of taxa: %i' % len(taxa))
sys.stdout.write('type\tortholog\t' + taxa[0] + '\t' + taxa[1] + '\n')
k = 0
for comps, cmptype in itertools.izip([o_o, o_m, m_m], ['1-1', '1-M', 'M-M']):
    for c in comps:
        c_seqs = {t:[] for t in taxa}
        for s in c:
            t = s.split('::')[0]
            c_seqs[t].append(s)
        for p in itertools.product(c_seqs[taxa[0]], c_seqs[taxa[1]]):
            sys.stdout.write(cmptype + '\t' + str(k) + '\t' + p[0] + '\t' + p[1] + '\n')
        k += 1
