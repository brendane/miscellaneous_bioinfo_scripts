#!/usr/bin/env python3
"""
    Process the long output from Ranger-DTL.

    ranger_dtl_parse_long_output.py --output <output file> <input file>
"""

import argparse
import collections
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

## Parse the long_output file
in_sp_t = False
in_g_t = False
in_reconcile = False
sp_tree = None
gene_trees = []
reconciliations = []
with open(args.input, 'rt') as ih:
    for line in ih:
        line = line.strip()
        if in_sp_t and line != '':
            sp_tree = line
            in_sp_t = False
        elif line == 'Species Tree:':
            in_sp_t = True
        elif in_g_t and line != '':
            gene_trees.append(line)
        elif in_g_t and line == '':
            in_g_t = False
        elif 'Optimal Rootings for Gene Tree' in line:
            in_g_t = True
        elif line.startswith('Reconciliation Data'):
            in_reconcile = True
        elif in_reconcile and line != '' and not line.startswith('Reconciliation Data'):
            reconciliations.append(line)

## Process the reconcilations
recon_table = collections.defaultdict(lambda: {'D':0, 'T':0, 'C':0})
for recon in reconciliations:
    sp_node = re.sub(',.+', '', re.sub('.+Most Frequent mapping --> ', '', recon))
    d = int(re.sub('.+Duplications = ([0-9]+),.+', '\\1', recon))
    s = int(re.sub('.+Speciations = ([0-9]+),.+', '\\1', recon))
    t = int(re.sub('.+Transfers = ([0-9]+).+', '\\1', recon))
    recon_table[sp_node]['C'] += (d + s + t)
    recon_table[sp_node]['D'] += d
    recon_table[sp_node]['T'] += t

with open(args.output, 'wt') as oh:
    oh.write('#' + sp_tree + '\n')
    oh.write('node\tduplications\ttransfers\n')
    for node, data in recon_table.items():
        oh.write(node + '\t' + str(data['D']/float(data['C'])) + '\t' +
                 str(data['T']/float(data['C'])) + '\n')
