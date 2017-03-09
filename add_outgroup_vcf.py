#!/usr/bin/env python
"""
    A simple script to add outgroup information to the INFO field of a
    VCF file.

    add_outgroup_vcf.py <outgroup file> <input VCF> <field ID>
"""

import gzip
import sys

og_file = sys.argv[1]
vcf_file = sys.argv[2]
id_name = sys.argv[3]

og_alleles = {}
with open(og_file, 'rb') as handle:
    for line in handle:
        row = line.strip().split('\t')
        og_alleles[(row[0], row[1])] = row[2].upper()

out = sys.stdout

open_fun = open
if vcf_file.endswith('.gz'):
    open_fun = gzip.open

data = False
with open_fun(vcf_file, 'rb') as handle:
    for line in handle:
        if line.startswith('#CHROM'):
            out.write('##INFO=<ID=' + id_name + ',Number=1,Type=Integer,' +
                      'Description="Whether reference matches outgroup">\n')
            out.write(line)
            data = True
        elif data:
            fields = line.strip().split('\t')
            ref_allele = fields[3].upper()
            alt_allele_len = max(len(x) for x in fields[4].split(','))
            chrom = fields[0]
            pos = fields[1]
            key = (fields[0], fields[1])
            if key in og_alleles and len(ref_allele) == 1 and alt_allele_len == 1:
                if og_alleles[key] == ref_allele:
                    ref_match_string = id_name + '=1'
                else:
                    ref_match_string = id_name + '=0'
            else:
                ref_match_string = id_name + '=-1'
            fields = line.strip().split('\t')
            info_field = fields[7]
            if info_field == '.':
                info_field = ref_match_string
            else:
               info_field += ';' + ref_match_string
            fields[7] = info_field
            out.write('\t'.join(fields) + '\n')
        else:
            out.write(line)
