#!/usr/bin/env python2
"""
    Get allele frequencies from PoPoolation synchronized file format.

    allele_freq_sync.py --output <output file> [--subset <subset>] 
        <pool name> <pool file> <sync files...>

    The pool file should list the names of the pools in order, with one
    pool per line.

    Allele frequencies calculations exclude Ns. The minor allele frequency
    is calculated within just the target pool, and it is calculated as
    1 - freq(most common allele). If there are more than two alleles,
    this could be > 0.5.

    Note also that these allele frequencies are just from the base counts
    - no fancy math here.
"""

#==============================================================================#

import argparse
import csv
import itertools

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--subset')
parser.add_argument('pool')
parser.add_argument('pfile')
parser.add_argument('sync', nargs='+')
args = parser.parse_args()

pool = args.pool

# If a subset of positions is desired
if args.subset is None:
    positions = None
else:
    positions = set()
    with open(args.subset, 'rb') as handle:
        for line in handle:
            contig, _, pos = line.strip().split()[:3]
            positions.add((contig, pos))

# Make a list of the pools in the sync file and get the index
# of the target pool.
pool_names = []
with open(args.pfile, 'rb') as inhandle:
    for line in inhandle:
        pool_names.append(line.strip())
pool_idx = pool_names.index(pool)

# Scan through the sync files and report the allele frequency
with open(args.output, 'wb') as out:
    out.write('contig\tpos\tref_allele\tref_freq\tmaf\tn_alleles\tn_reads\n')
    for sync_file in args.sync:
        with open(sync_file, 'rb') as inhandle:
            rdr = csv.reader(inhandle, delimiter='\t')
            for row in rdr:
                contig = row[0]
                pos = row[1]
                if (positions is not None) and ((contig, pos) not in positions):
                    continue
                ref_allele = row[2].upper()
                counts = {}
                for allele, c in itertools.izip(['A','T','C','G','N','D'],
                                                row[pool_idx + 3].split(':')):
                    if allele == 'N':
                        continue
                    counts[allele] = int(c)
                total = float(sum(counts.itervalues()))
                if total == 0:
                    ref_freq = 0
                    minor_freq = 0.0
                    n_alleles = 0
                else:
                    ref_freq = counts[ref_allele] / total
                    major_count = max(counts.itervalues())
                    minor_freq = (total - major_count) / total
                    n_alleles = sum(1 for c in counts.itervalues() if c > 0)
                out.write(contig + '\t' + pos + '\t' + ref_allele + '\t' +
                          str(ref_freq) + '\t' + str(minor_freq) + '\t' +
                          str(n_alleles) + '\t' + str(total) + '\n')

