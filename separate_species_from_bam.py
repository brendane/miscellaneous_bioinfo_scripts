#!/usr/bin/env python2.7
"""
    Given a bam file and a list of contigs/chrs, separate reads that
    align to one set of contigs/chrs from the other reads.

    separate_species_from_bam.py --output <output file prefix>
        --first <first species name> --second <second species name>
        <bam file> <comma-separated list of contigs from the first species>

    Input bam file should be sorted by read name or there will be
    complications with the output.

    Output files are interleaved.
"""

import argparse

import pysam as ps

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--first')
parser.add_argument('--second')
parser.add_argument('bam')
parser.add_argument('contigs')
args = parser.parse_args()

outpre = args.output
spp0 = args.first
spp1 = args.second
if spp0 is None or spp1 is None:
    parser.print_help()
    exit(1)
bam_file = args.bam
contigs = set(args.contigs.split(','))

# Figure out which species each read aligns to
# Note that I'm assuming that qname does not include the /1 or /2
# part of the name (F or R read), so that this will count mates as
# a single read.
aln_0 = set()
aln_1 = set()
unmapped_count = 0
for read in ps.Samfile(bam_file, 'rb'):
    if read.is_unmapped:
        unmapped_count += 1
        continue
    if read.reference_name in contigs:
        aln_0.add(read.qname)
    else:
        aln_1.add(read.qname)
only_0 = aln_0 - aln_1
only_1 = aln_1 - aln_0
all_reads = set.union(aln_0, aln_1)
print '%i reads aligned to %s'      % (len(aln_0), spp0)
print '%i reads aligned to %s'      % (len(aln_1), spp1)
print '%i reads not aligned'        % (unmapped_count)
print '%i reads aligned in total'   % (len(all_reads))
print '%i reads aligned only to %s (%f)' % (len(only_0), spp0, len(only_0) / round(float(len(aln_0)), 2))
print '%i reads aligned only to %s (%f)' % (len(only_1), spp1, len(only_1) / round(float(len(aln_1)), 2))

# Write reads that align to only one species to fastq files
# Both mates in a pair are either written to one file or not written
# to anything.
written = set()
with open(outpre + '.' + spp0 + '.fastq', 'wb') as out0:
    with open(outpre + '.' + spp1 + '.fastq', 'wb') as out1:
        for read in ps.Samfile(bam_file, 'rb'):
            qn = read.qname
            if read.is_unmapped and qn not in only_0 and qn not in only_1:
                continue
            name = read.qname + '/' + ('1' if read.is_read1 else '2')
            if name not in written and qn in only_0:
                out0.write('@' + name + '\n' + read.seq + '\n+\n' +
                           read.qual + '\n')
                written.add(name)
            elif name not in written and qn in only_1:
                out1.write('@' + name + '\n' + read.seq + '\n+\n' +
                           read.qual + '\n')
                written.add(name)
            else:
                pass
