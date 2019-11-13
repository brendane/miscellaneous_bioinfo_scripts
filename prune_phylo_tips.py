#! /usr/bin/env python

import sys

import dendropy

tree_str = ''
with open(sys.argv[1], 'rt') as ih:
    for line in ih:
        tree_str += line.strip()

include = []
with open(sys.argv[2], 'rt') as ih:
    for line in ih:
        include.append(line.strip())

tree = dendropy.Tree.get(data=tree_str, schema="newick")
tree.retain_taxa_with_labels(include)
print(tree.as_string(schema='newick'))
