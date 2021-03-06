#!/usr/bin/env python3
"""
    Calculate betweenness centrality fo plant and bacterial genes using
    networkx. Relies on having a very specific input data structure.

    betweenness_centrality_plant_bacteria.py <adjacency file> <threshold>
"""

import csv
import sys

import networkx as nx
from networkx.algorithms.centrality.betweenness_subset import betweenness_centrality_subset as bcs
from networkx.algorithms.centrality.betweenness import betweenness_centrality as bc

thresh = float(sys.argv[2])

## Make the graph
g = nx.Graph()
with open(sys.argv[1], 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for i, row in enumerate(rdr):
        if i == 0:
            genes = row
            plant_genes = [x for x in row if x.startswith('Med')]
            bact_genes = [x for x in row if not x.startswith('Med')]
            for gene in genes:
                g.add_node(gene)
        else:
            focal_gene = row[0]
            weights = [float(x) for x in row[1:]]
            for j, (gene, weight) in enumerate(zip(genes, weights)):
                if j < i:
                    ## Note: first row of file is the header, so j < i, not j <= i
                    continue
                if weight < thresh:
                    continue
                d = 1-weight
                g.add_edge(focal_gene, gene, weight=weight, dist=d)
sys.stderr.write('Graph construction finished\n')

sys.stdout.write('gene\tcomparison\tbetweenness\n')
bb = bcs(g, sources=bact_genes, targets=bact_genes,
        normalized=False, weight='dist')
for gene in genes:
    sys.stdout.write(gene + '\tbacteria\t' + str(bb[gene]) + '\n')
bp = bcs(g, sources=plant_genes, targets=plant_genes,
        normalized=False, weight='dist')
for gene in genes:
    sys.stdout.write(gene + '\tplant\t' + str(bp[gene]) + '\n')
b = bc(g, normalized=False, weight='dist')
for gene in genes:
    sys.stdout.write(gene + '\tall\t' + str(b[gene]) + '\n')
