#!/usr/bin/env python2
"""
    Choose one variant to represent each group from a file with LD
    groupings. Chooses the first one in the file, so the way the file
    is sorted will affect the choice.

    ld_grouping_choose_one.py --output <output file> <input file>

    The output is in bed format, with columns five and six holding the
    group number and all the variants associated with the group,
    respectively.
"""

#==============================================================================#

import argparse
import collections
import csv

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

choices = {}
groupings = collections.defaultdict(list)
with open(args.input, 'rb') as inhandle:
    rdr = csv.DictReader(inhandle, delimiter='\t')
    for row in rdr:
        key = (row['chrom'], int(row['pos']))
        grp = int(row['group'])
        if grp not in choices:
            choices[grp] = key
        groupings[grp].append(row)

with open(args.output, 'wb') as out:
    for g, k in sorted(choices.iteritems()):
        out.write(k[0] + '\t' + str(k[1]-1) + '\t' + str(k[1]) + '\t' +
                  k[0] + '-' + str(k[1]) + '\t' + str(g) + '\t')
        out.write(','.join(x['chrom'] + '-' + str(x['pos'])
                           for x in groupings[g]) + '\n')

