#!/usr/bin/env python3
"""
    Given a rooted phylogenetic tree in newick format, find sets of taxa
    (tips) that are monophyletic and no bigger than a certain number
    of taxa. For closely related clusters, will split haphazardly,
    ignoring phylogenetic relationships. This makes the chunks
    more evenly-sized.

    break_tree_into_chunks.py --output <output prefix>
        <tree> <max chunk size> <max tree distance to split node>
"""

import argparse
import itertools
import math

import dendropy

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    for x in itertools.zip_longest(*args, fillvalue=fillvalue):
        yield [xx for xx in x if xx is not None]

def get_leaf_taxa(node, targets=None):
    if node.is_leaf():
        genes = {node.taxon.label}
    else:
        genes = {n.taxon.label for n in node.leaf_nodes()}
    if targets:
        genes = genes.intersection(targets)
    return genes

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('tree')
parser.add_argument('chunk', type=int)
parser.add_argument('scalar', type=float)
args = parser.parse_args()

tree = dendropy.Tree.get(path=args.tree, schema='newick',
                         rooting='default-rooted',
                         preserve_underscores=True)

## Get distances
distance_matrix = tree.phylogenetic_distance_matrix()
max_distance = 0
for t1 in tree.taxon_namespace:
    for t2 in tree.taxon_namespace:
        d = distance_matrix(t1, t2)
        if d > max_distance:
            max_distance = d

combine_threshold = d * args.scalar

## Find the largest chunks of the tree that are monophyletic and
## are no more than about args.chunk in size
chunks = []
for node in tree.preorder_node_iter():
    tips = get_leaf_taxa(node)
    skip = False
    if len(chunks) > 0:
        for c in chunks:
            if len(set.intersection(c, tips)) > 0:
                ## Child of already added split, continue
                skip = True
                break
    if skip: continue
    if len(tips) <= args.chunk:
        chunks.append(tips)
    elif len(tips) <= args.chunk * 10:
        md = 0
        for t1 in tree.taxon_namespace.findall(list(tips)):
            for t2 in tree.taxon_namespace.findall(list(tips)):
                d = distance_matrix(t1, t2)
                if d > md:
                    md = d
        if md <= combine_threshold:
            ## All of the tips under this node are close to
            ## each other, so just split the node into a few
            ## large pieces, while ignoring phylogenetic relationships
            for tp in grouper(tips, math.ceil(len(tips) / math.ceil(len(tips) / args.chunk))):
                chunks.append(set(tp))

## We don't want chunks with only 1 strain, so combine based on
## phylogenetic similarity. There may still be a taxon left over
## at the end.
single_chunks = [chunk for chunk in chunks if len(chunk) == 1]
combinations = []
done = set()
for sc in single_chunks:
    t = list(sc)[0]
    if t in done:
        continue
    closest_distance = None
    closest = None
    tn = tree.taxon_namespace.findall(t)[0]
    for sc2 in single_chunks:
        t2 = list(sc2)[0]
        if t == t2:
            continue
        if t2 in done:
            continue
        tn2 = tree.taxon_namespace.findall(t2)[0]
        d = distance_matrix(tn, tn2)
        if closest_distance is None or d < closest_distance:
            closest_distance = d
            closest = t2
    if closest is not None:
        combinations.append({t, closest})
        done.add(t); done.add(closest)
    else:
        combinations.append({t})
        done.add(t)

final_chunks = [chunk for chunk in chunks if len(chunk) > 1]
final_chunks += combinations

## Write files ending in .*.txt where * = chunk number (from zero)
for i, chunk in enumerate(final_chunks):
    with open(args.output + '.' + str(i) + '.txt', 'wt') as oh:
        for taxon in chunk:
            oh.write(taxon + '\n')
