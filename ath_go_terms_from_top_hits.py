#!/usr/bin/env python3
"""
    Take tab-delimited file and table of GO terms from
    arabidopsis.org to annotate top hits with GO terms.

    ath_go_terms_from_top_hits.py <top hits> <go terms>

    Input file with top hits should have:
    1. query
    2. gene in Arabidopsis
    3. bitscore
    4. evalue
"""

import collections
import csv
import gzip
import re
import sys

go = collections.defaultdict(set)
all_go_terms = set()
open_fun = open
if sys.argv[2].endswith('.gz'): open_fun = gzip.open
with open_fun(sys.argv[2], 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        gene = row[0]
        term_number = row[5]
        description = row[3] + ' ' + row[4]
        go[gene].add((term_number, description))
        all_go_terms.add((term_number, description))

sys.stdout.write('query\tlocus\tscore\tevalue\tGO\tGO_description\n')
with open(sys.argv[1], 'r') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        gene = re.sub('\\.[0-9]+$', '', row['hit'])
        go_terms = go[gene]
        for g in go_terms:
            sys.stdout.write(row['query'] + '\t' +
                             gene + '\t' +
                             row['score'] + '\t' +
                             row['evalue'] + '\t' +
                             g[0] + '\t' + g[1] + '\n')
