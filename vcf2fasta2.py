#!/usr/bin/env python2
"""
    Using a VCF file, a reference genome, and a table of read-depths,
    make a FASTA file for a list of individuals.

    vcf2fasta.pyt --output <output file> [OPTIONS] <vcf> <reference>
        <depth file>

    --min-coverage      Minimum depth to count a site as covered
    --replicon          Replicon to include
    --maf               Minimum minor allele frequency

    Assumes only bi-allelic SNPs in VCF file, and only works on haploids.

    UPDATE 05 March 2018: Can now work on a single replicon.

    UPDATE 08 June 2018: Added minor allele frequency filter
"""

#==============================================================================#

import argparse
import copy
import csv

from Bio import SeqIO
import pysam

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--replicon')
parser.add_argument('--min-coverage', type=int)
parser.add_argument('--max-pos', type=int, default=1E9)
parser.add_argument('--maf', type=float, default=0.)
parser.add_argument('vcf')
parser.add_argument('ref')
parser.add_argument('depth')
args = parser.parse_args()

vcf_file_name = args.vcf
ref_file_name = args.ref
depth_file_name = args.depth
output_file_name = args.output
min_coverage = args.min_coverage
max_pos = args.max_pos
min_maf = args.maf
replicon_to_use = args.replicon

# Read reference sequence
ref_seq = {}
for rec in SeqIO.parse(ref_file_name, 'fasta'):
    if replicon_to_use is not None and rec.id != replicon_to_use:
        continue
    ref_seq[rec.id] = list(str(rec.seq))

# Get variant positions from VCF file
vf = pysam.VariantFile(vcf_file_name)
strains = list(vf.header.samples)
snps = {strain:{} for strain in strains}
for rec in vf:
    # Skip indels
    if max(map(len, rec.alleles)) > 1:
        continue
    ref_pos = rec.pos
    ref_replicon = rec.chrom
    if ref_replicon not in ref_seq:
        continue
    if replicon_to_use is not None and ref_replicon != replicon_to_use:
        continue
    if ref_pos > max_pos:
        continue
    ref_allele = None
    ref_count = 0
    alt_count = 0
    for sample in rec.samples:
        sample_data = rec.samples[sample]
        gt = sample_data.alleles[0]
        if gt is None or (len(sample_data['GT']) > 1 and sample_data.alleles[1] != gt):
            gt = 'N'
        if gt != 'N':
            if ref_allele is None or gt == ref_allele:
                ref_allele = gt
                ref_count += 1
            else:
                alt_count += 1
        snps[sample][(ref_replicon, ref_pos)] = gt
    if min_maf > 0.:
        maf = float(ref_count) / (alt_count + ref_count)
        if maf < min_maf or (1 - maf) < min_maf:
            for sample in rec.samples:
                snps[sample][(ref_replicon, ref_pos)] = 'N'
        

# Process the depth file
not_covered = {s:set() for s in strains}
with open(depth_file_name, 'rb') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    ninds = len(rdr.fieldnames) - 2
    for row in rdr:
        for strain in strains:
            covered = row[strain] > min_coverage
            if not covered:
                not_covered[strain].add((row['replicon'], int(row['pos'])))

# Now write the fasta file
with open(output_file_name, 'wb') as out:
    for strain in strains:
        out.write('>' + strain + '\n')
        sequence = copy.deepcopy(ref_seq)
        for r, p in not_covered[strain]:
            sequence[r][p-1] = 'N'
        for (r, p), gt in snps[strain].iteritems():
            sequence[r][p-1] = gt
        seq = []
        for replicon_sequence in sequence.itervalues():
            seq += replicon_sequence
        out.write(''.join(seq) + '\n')
