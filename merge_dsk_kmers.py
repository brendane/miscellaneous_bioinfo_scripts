#!/usr/bin/env python2.7
"""
    Merge kmer counts from dsk.

    merge_dsk_kmers.py [--count] [--factor] --output <output prefix>
        <input file list>
"""

import argparse
import collections
import gzip
import os

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--count', type=int, default=-1)
parser.add_argument('--factor', type=float, default=-1)
parser.add_argument('--output')
parser.add_argument('input')
parser.add_argument('min_strain_count', type=int)
args = parser.parse_args()

min_strain_count = args.min_strain_count
min_count = None
min_factor = None
if args.count >= 0:
    min_count = args.count
elif args.factor >= 0:
    min_factor = args.factor
else:
    raise Exception('Have to provide minimum count or factor')

outpre = args.output

# Get the list of files to read and strain names
infiles = {}
with open(args.input, 'rb') as ih:
    for line in ih:
        s, f = line.strip().split('\t')
        infiles[f] = s

# Get the list of all kmers
kmers = collections.defaultdict(int)
means = {}
for infile in infiles:
    of = open
    if infile.endswith('.gz'):
        of = gzip.open
    with of(infile, 'rb') as ih:
        c = 0.
        s = 0
        for line in ih:
            kmers[(line.strip().split()[0])] += 1
            s += int(line.strip().split()[1])
            c += 1.
    means[infiles[infile]] = s / c

remove = set()
for k, c in kmers.iteritems():
    if c < min_strain_count:
        remove.add(k)
for k in remove:
    del kmers[k]
del remove

# Then check presence or absence of each kmer in each strain, save
# to temporary files
tmp_files = {}
for infile, strain in infiles.iteritems():
    kmer_counts = {k:0 for k in kmers}
    of = open
    if infile.endswith('.gz'):
        of = gzip.open
    with of(infile, 'rb') as ih:
        for line in ih:
            k, c = line.strip().split()
            kmer_counts[k] = int(c)
    tmp = outpre + '.' + strain
    tmp_files[strain] = tmp
    with open(tmp, 'wb') as oh:
        for k in kmers:
            if min_count is None:
                present = kmer_counts[k] >= (means[strain] * min_factor)
            else:
                present = kmer_count[k] >= min_count
            oh.write(k + '\t' + str(int(present)) + '\n')

"""
# Put everything into one file; has to open lots of file handles
handles = {s:open(f, 'rb') for s, f in tmp_files.iteritems()}
with open(outfile, 'wb') as oh:
    oh.write('kmer')
    for s in handles:
        oh.write('\t' + s)
    oh.write('\n')
    for k in kmers:
        x = []
        c = 0
        for s in handles:
            kk, pa = handles[s].readline().strip().split('\t')
            if kk != k:
                raise Exception('Kmer doesn\'t match')
            x.append(pa)
            if pa == '1':
                c += 1
        if c > 0:
            oh.write(k + '\t' + '\t'.join(x) + '\n')

# Close handles and delete temporary files
for h in handles.itervalues():
    h.close()
for f in tmp_files.itervalues():
    os.unlink(f)
"""
