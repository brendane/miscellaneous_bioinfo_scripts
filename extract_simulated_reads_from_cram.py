#!/usr/bin/env python3
"""
    Given simulated strain frequencies and reads from individual
    sequencing, create a cram file with sampled reads.

    extract_simulated_reads_cram.py --output <output prefix>
            --nreads <number of total reads (approx.)>
            --slop <error in number of reads>
            <file with read locations> <file with frequencies>

    Requires biopython and pysam. Assumes paired end reads.
"""

import argparse
import gzip
import math
import random

from Bio import SeqIO
import pysam

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--nreads', type=int)
parser.add_argument('--slop', type=float, default=0.)
parser.add_argument('aligns')
parser.add_argument('freqs')
args = parser.parse_args()

## Read in the file locations and format for each strain
## This file should be tab-separated, no header, with columns:
## strain, file
print('reading file locations')
files = {}
with open(args.aligns, 'r') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        files[fields[0]] = fields[1]

## Read in frequencies. This file should be tab-separated, no header,
## with columns: strain, frequency.
print('reading frequencies file')
freqs = {}
with open(args.freqs, 'r') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        freqs[fields[0]] = float(fields[1])

## Calculate the number of reads to take from each file
print('calculating number of reads')
rand = random.Random()
n_reads = {}
for strain, freq in freqs.items():
    slop = rand.uniform(1-args.slop, 1+args.slop)
    n_reads[strain] = int(math.ceil(freq * args.nreads * slop))

## Get the reads
print('sampling reads')
first = True
rand = random.Random()
oh = None
nr = {}
for strain, fname in files.items():
    if strain not in n_reads:
        continue
    nr[strain] = 0
    with pysam.AlignmentFile(fname, 'rc') as ih:
        reads_sampled = set()
        # 1. If first, get header and apply to output file
        if first:
            oh = pysam.AlignmentFile(args.output + '.cram', 'wc', template=ih)
            first = False
        # 2. Calculate number of reads in file total; divide by two to get pairs
        total_reads = ih.count() / 2
        prob = total_reads / float(n_reads[strain])
        # 3. Randomly sample with appropriate probability
        for read in ih.fetch():
            if read.qname in reads_sampled:
                # Write this read -- already sampled mate
                oh.write(read)
                pass
            elif rand.uniform(0, 1) < prob:
                # Passes probability filter
                reads_sampled.add(read.qname)
                oh.write(read)
                nr[strain] += 1
oh.close()

with open(args.output + '.nreads.txt', 'w') as nr:
    for strain, n in nr.items():
        nr.writelines([strain, '\t', n, '\n'])
