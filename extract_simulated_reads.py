#!/usr/bin/env python2.7
"""
    Given simulated strain frequencies and reads from individual
    sequencing, create a pair of fastq files with the correct number
    of reads.

    extract_simulated_reads.py --output <output prefix>
            --nreads <number of total reads (approx.)>
            --length <length of reads>
            --slop <error in number of reads>
            <file with read locations> <file with frequencies>

    A random number between 1-slop and 1+slop is multiplied by the
    number of reads that are supposed to be sampled from each file
    to add some error. Set to zero for no error.

    Requires biopython.
"""

import argparse
import gzip
import itertools
import math
import random

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--nreads', type=int)
parser.add_argument('--length', type=int, default=0)
parser.add_argument('--slop', type=float, default=0.)
parser.add_argument('reads')
parser.add_argument('freqs')
args = parser.parse_args()

## Read in the file locations and format for each strain
## This file should be tab-separated, no header, with columns:
## strain, format, reads 1 file, reads 2 file
files = {}
with open(args.reads, 'rb') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        files[fields[0]] = tuple(fields[1:])

## Read in frequencies. This file should be tab-separated, no header,
## with columns: strain, frequency.
freqs = {}
with open(args.freqs, 'rb') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        freqs[fields[0]] = float(fields[1])

## Calculate the number of reads to take from each file
rand = random.Random()
n_reads = {}
for strain, freq in freqs.iteritems():
    slop = rand.uniform(1-args.slop, 1+args.slop)
    n_reads[strain] = int(math.ceil(freq * args.nreads * slop))

## Get the reads
with open(args.output + '.r1.fastq', 'wb') as o1:
    with open(args.output + '.r2.fastq', 'wb') as o2:
        for strain, file_info in files.iteritems():
            if strain not in n_reads:
                continue
            ofun = open
            if file_info[1].endswith('.gz'):
                ofun = gzip.open
            rnames = []
            with ofun(file_info[1]) as r1h:
                for rec in SeqIO.parse(r1h, file_info[0]):
                    rnames.append(rec.id)
            chosen = set()
            while len(chosen) < n_reads[strain]:
                chosen.add(rand.choice(rnames))
            with ofun(file_info[1]) as r1h:
                with ofun(file_info[2]) as r2h:
                    p1 = SeqIO.parse(r1h, file_info[0])
                    p2 = SeqIO.parse(r2h, file_info[0])
                    for r1, r2 in itertools.izip(p1, p2):
                        if r1.id in chosen:
                            if args.length > 0:
                                r1 = r1[:args.length]
                                r2 = r2[:args.length]
                            SeqIO.write(r1, o1, 'fastq')
                            SeqIO.write(r2, o2, 'fastq')
