#!/usr/bin/env python2
"""
    Take the output of xmfa2outgroup.py and create a fasta file with
    the outgroup sequence. Missing or ambiguous sites are coded as "N."

    outgroup2fasta.py --output <output file> --replicon <target replicon>
        <input> <fai file> <sequence name>
"""

#==============================================================================#

import argparse
import csv

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--replicon')
parser.add_argument('input')
parser.add_argument('fai')
parser.add_argument('name')
args = parser.parse_args()

replicon = args.replicon

replicon_lengths = {}
with open(args.fai, 'rb') as handle:
    for line in handle:
        fields = line.strip().split()
        replicon_lengths[fields[0]] = int(fields[1])

bases = ['N'] * replicon_lengths[replicon]
with open(args.input, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        if row[0] != replicon:
            continue
        bases[int(row[1])-1] = row[2].upper()

with open(args.output, 'wb') as out:
    out.write('>' + args.name + '\n')
    out.write(''.join(bases) + '\n')

