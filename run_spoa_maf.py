#!/usr/bin/env python3
"""
    Use spoa to get a MAF file chunk.

    run_spoa_maf.py [SPOA OPTIONS] <gff dir> <fasta dir> <fai> <label list>
"""

import argparse
import re
import subprocess
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('-l', default='1')
parser.add_argument('-r', default='1')
parser.add_argument('-e', default='-8')
parser.add_argument('gff')
parser.add_argument('fasta')
parser.add_argument('fai')
parser.add_argument('label')
args = parser.parse_args()

## TODO: modify to loop through a gff file given a list of IDs

## Source lengths
lens = {}
with open(args.fai, 'rt') as ih:
    for line in ih:
        row = line.strip().split('\t')
        lens[row[0]] = int(row[1])

labels = []
with open(args.label, 'rt') as ih:
    for line in ih:
        labels.append(line.strip())

for label in labels:
    ## Make alignment
    sys.stderr.write('LCB=' + label + '\n')
    p = subprocess.Popen(['spoa', '-r', args.r, '-l', args.l,
                          '-e', args.e,
                          args.fasta + '/' + label + '.fasta'],
                         stdout=subprocess.PIPE)
    aln = list(map(lambda x: re.sub(':.+', '', x),
                   p.communicate()[0].decode('utf-8').split('\n')[1:]))

    ## Get sequence lengths
    lengths = []
    for s in aln:
        lengths.append(len(s) - s.count('-'))

    ## Get sequence directions, IDs, and start/end
    strands = []
    ids = []
    starts = []
    ends = []
    with open(args.gff + '/' + label + '.gff', 'rt') as ih:
        for line in ih:
            if line.startswith('#'): continue
            row = line.strip().split('\t')
            lab = re.sub('id=', '', row[8])
            strands.append(row[6])
            ids.append(row[0])
            starts.append(int(row[3])) # 1-based
            ends.append(int(row[4]))

    ## Get src lengths
    src_lengths = []
    for i in ids:
        src_lengths.append(lens[i])

    ## Put together
    sys.stdout.write('a label=%s\n' % lab)
    for i, s, e, st, sl, l, a in zip(ids, starts, ends, strands, src_lengths, lengths, aln):
        if st == '+':
            sys.stdout.write('\t'.join(map(str,
                                           ['s', i, s-1, l, st, sl, a])))
            sys.stdout.write('\n')
        else:
            start = sl - e
            sys.stdout.write('\t'.join(map(str,
                                           ['s', i, start, l, st, sl, a])))
            sys.stdout.write('\n')
    sys.stdout.write('\n')
