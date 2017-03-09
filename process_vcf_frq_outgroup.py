#!/usr/bin/env python
"""
    A simple script to add modify the frq output from VCFtools and add
    outgroup information.

    process_vcf_frq_outgroup.py <outgroup file> <input .frq file>
"""

import gzip
import sys

og_file = sys.argv[1]
frq_file = sys.argv[2]

og_alleles = {}
with open(og_file, 'rb') as handle:
    for line in handle:
        row = line.strip().split('\t')
        og_alleles[(row[0], row[1])] = row[2].upper()

out = sys.stdout

out.write('chrom\tpos\tn_alleles\tn_strains\tref_allele\tref_freq\t' +
          'outgroup_allele\toutgroup_match\n')
with open(frq_file, 'rb') as handle:
    for line in handle:
        if line.startswith('CHROM'):
            continue
        fields = line.strip().split()
        chrom, pos, n_all, n_chr, ref_af_str = fields[:5]
        ref_all, ref_freq = ref_af_str.split(':')
        ref_all = ref_all.upper()
        if len(ref_all) == 1 and (chrom, pos) in og_alleles:
            og_all = og_alleles[(chrom, pos)]
            if og_all == ref_all:
                match = '1'
            else:
                match = '0'
        else:
            og_all = 'N'
            match = 'NA'
        out.write(chrom + '\t' + pos + '\t' + n_all + '\t' +
                  str(int(n_chr)/2) + '\t' + ref_all + '\t' +
                  ref_freq + '\t' + og_all + '\t' + match + '\n')
