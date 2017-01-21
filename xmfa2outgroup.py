#!/usr/bin/env python
"""
    Take an xmfa file from mugsy with a pairwise alignment
    and generate outgroup states.

    xmfa2outgroup.py --output <output> <input> <outgroup>
"""

#==============================================================================#

import argparse
import itertools

from Bio import SeqIO

#==============================================================================#

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
parser.add_argument('outgroup')
args = parser.parse_args()

outgroup = args.outgroup
outgroup_states = {}

# strain.replicon 0-based-start align-len strand replicon-len
#
# if strand == '-':
#   The start coordinate is still the beginning of the sequence,
#   but the sequence is complemented

for rec1, rec2 in grouper(SeqIO.parse(args.input, 'fasta'), 2):
    if rec2[-1] != '=':
        raise Exception('Alignment is not pairwise')
    rec2 = rec2[:-1] # get rid of "="
    strain1 = rec1.description.split(' ')[0].split('.')[0]
    strain2 = rec2.description.split(' ')[0].split('.')[0]
    if strain1 == outgroup:
        o = rec1
        i = rec2
    elif strain2 == outgroup:
        o = rec2
        i = rec1
    else:
        raise Exception('outgroup strain not found')

    o_start = int(o.description.split(' ')[1])
    o_strand = o.description.split(' ')[3]
    i_start = int(i.description.split(' ')[1])
    i_strand = i.description.split(' ')[3]
    i_repl = i.description.split(' ')[0].split('.')[1]

    j = 0
    for ib, ob in itertools.izip(i, o):
        if ib == '-':
            continue
        p = j + i_start + 1 # convert to 1-based
        if i_strand == '-':
            ob = Seq.reverse_complement(ob).upper()
        else:
            ob = ob.upper()
        key = (i_repl, p)
        if key in outgroup_states and outgroup_states[key] != ob:
            outgroup_states[key] = 'N'
        elif ob != '-':
            outgroup_states[key] = ob
        j += 1

with open(args.output, 'wb') as out:
    for (repl, p), o in sorted(outgroup_states.iteritems()):
        if o == 'N':
            continue
        out.write(repl + '\t' + str(p) + '\t' + o + '\n')
