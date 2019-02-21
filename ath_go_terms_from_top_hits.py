#!/usr/bin/env python3
"""
    Take tab-delimited file and table of GO terms from
    arabidopsis.org to annotate top hits with GO terms.

    ath_go_terms_from_top_hits.py <top hits> <go terms> [<go term category>]

    Input file with top hits should have:
    1. query
    2. gene in Arabidopsis
    3. bitscore
    4. evalue

    Reports GO Slim terms.
"""

import argparse
import collections
import csv
import gzip
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--aspect')
parser.add_argument('--slim', default=False, action='store_true')
parser.add_argument('hits')
parser.add_argument('terms')
args = parser.parse_args()

cat = args.aspect

go = collections.defaultdict(set)
open_fun = open
with open(args.terms, 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if cat is not None and cat != row[7]:
            continue
        gene = row[0]
        term_number = row[5]
        #description = row[3] + ' ' + row[4]
        if args.slim:
            description = row[8]
            go[gene].add(description)
        else:
            description = row[4]
            go[gene].add((term_number, description))

if args.slim:
    sys.stdout.write('query\tlocus\tscore\tevalue\tGO_description\n')
else:
    sys.stdout.write('query\tlocus\tscore\tevalue\tGO\tGO_description\n')
with open(args.hits, 'r') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        gene = re.sub('\\.[0-9]+$', '', row['hit'])
        go_terms = go[gene]
        for g in go_terms:
            if args.slim:
                gg = g
            else:
                gg = g[0] + '\t' + g[1]
            sys.stdout.write(row['query'] + '\t' +
                             gene + '\t' +
                             row['score'] + '\t' +
                             row['evalue'] + '\t' +
                             gg + '\n')
