#!/usr/bin/env python2
"""
    Using a VCF file, a reference genome, and a table of read-depths,
    make a FASTA file for an individual.

    vcf2fasta.pyt --output <output file> [OPTIONS] <vcf> <reference>
        <depth file>

    --min-coverage      Minimum depth to count a site as covered
    --min-covered       Minimum proportion of strains covered to
                        count a site
    --sample            Which sample to target (this script does one at
                        a time)
    --replicon          Target replicon
    --het-missing       Set heterozygous calls to missing; otherwise, just
                        takes the first allele listed
    --max-pos           Maximum allowed position in replicon

    Assumes only bi-allelic SNPs in VCF file, and only works on haploids.
"""

#==============================================================================#

import argparse
import csv

from Bio import SeqIO
import pysam

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--min-coverage', type=int)
parser.add_argument('--min-covered', type=float)
parser.add_argument('--sample')
parser.add_argument('--replicon')
parser.add_argument('--het-missing', default=False, action='store_true')
parser.add_argument('--max-pos', type=int, default=1E9)
parser.add_argument('vcf')
parser.add_argument('ref')
parser.add_argument('depth')
args = parser.parse_args()

vcf_file_name = args.vcf
ref_file_name = args.ref
depth_file_name = args.depth
target_replicon = args.replicon
target_sample_name = args.sample
output_file_name = args.output
min_coverage = args.min_coverage
min_covered = args.min_covered
max_pos = args.max_pos

# Get variant positions from VCF file
snps = {}
vf = pysam.VariantFile(vcf_file_name)
vf.subset_samples([target_sample_name])
for rec in vf:
    ref_pos = rec.pos
    ref_replicon = rec.chrom
    if ref_replicon != target_replicon:
        continue
    if ref_pos > max_pos:
        continue
    ref_allele = rec.ref
    alt_allele = rec.alts[0]
    # Skip indels
    if max(map(len, [ref_allele] + list(rec.alts))) > 1:
        continue
    gt = (rec.samples[0]['GT'][0])
    if gt is None:
        snps[(ref_replicon, ref_pos)] = 'N'
    else:
        if args.het_missing and rec.samples[0]['GT'][1] != gt:
            snps[(ref_replicon, ref_pos)] = 'N'
        else:
            snps[(ref_replicon, ref_pos)] = ([ref_allele] + list(rec.alts))[gt]

# Read reference sequence
refseq = str(SeqIO.read(ref_file_name, 'fasta').seq)

# Process the depth file
not_covered = set()
with open(depth_file_name, 'rb') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    ninds = len(rdr.fieldnames) - 2
    for row in rdr:
        if row['replicon'] != target_replicon:
            continue
        covered = row[target_sample_name] > min_coverage
        if covered:
            n_covered = sum(1 for f, c in row.iteritems()
                            if f not in {'pos', 'replicon'} and
                            int(c) >= min_coverage)
            covered = covered and (float(n_covered) / ninds) >= min_covered
        if not covered:
            not_covered.add((row['replicon'], int(row['pos'])))

# Now write the fasta file
with open(output_file_name, 'wb') as out:
    out.write('>' + target_sample_name + '\n')
    for i in xrange(len(refseq)):
        key = (target_replicon, i)
        if key in not_covered:
            out.write('N')
        elif key in snps:
            out.write(snps[key].upper())
        else:
            out.write(refseq[i].upper())
    out.write('\n')
