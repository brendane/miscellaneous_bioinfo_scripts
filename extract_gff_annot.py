#!/usr/bin/env python3
"""
    Given a column in a tab-delimited file, extract certain gff
    formatted annotations. Prints the entire row with the extracted
    annotations substituted.

    extract_gff_annot.py <column number> <comma-separated fields> [<input>]

    column number is 1-based
"""

import sys
import urllib.parse

c = int(sys.argv[1])
fields = sys.argv[2].split(',')
if len(sys.argv) > 3:
    handle = open(sys.argv[3], 'rt')
else:
    handle = sys.stdin

with handle as ih:
    for line in ih:
        if line.startswith('#'):
            continue
        row = line.strip().split('\t')
        if '=' in row[c-1]:
            a = {x.split('=')[0]:x.split('=')[1] for x in row[c-1].split(';')}
        else:
            a = {}
        output = []
        for field in fields:
            if field in a:
                output.append(urllib.parse.unquote(a[field]))
            else:
                output.append('')
        row[c-1] = ' '.join(output)
        sys.stdout.write('\t'.join(row) + '\n')
