#!/usr/bin/env python2
"""
    From a tab-delimited file from snpEff, extract the most deleterious
    mutation type for each variant.

    Uses the mutation order from:
    http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf

    Changes upstream_gene_variant and downstream_gene_variant to
    intergenic and removes "_variant" from all the names.

    extract_worst_mutation.py --output <output file> <input file>

    If input or output is given as "-", read from / write to stdin/stdout.

    The input file should have the annotations in the third column.
"""

#==============================================================================#

import argparse
import re
import sys

#==============================================================================#

variant_order = [ 'chromosome_number_variation', 'exon_loss_variant',
                  'frameshift_variant', 'stop_gained', 'stop_lost',
                  'start_lost', 'splice_acceptor_variant', 'splice_donor_variant',
                  'rare_amino_acid_variant', 'missense_variant', 'inframe_insertion',
                  'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion',
                  '5_prime_UTR_truncation+exon_loss_variant', '3_prime_UTR_truncation+exon_loss',
                  'splice_branch_variant', 'splice_region_variant', 'splice_branch_variant',
                  'stop_retained_variant', 'initiator_codon_variant', 'synonymous_variant',
                  'initiator_codon_variant+non_canonical_start_codon',
                  'stop_retained_variant', 'coding_sequence_variant', '5_prime_UTR_variant',
                  '3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant',
                  'upstream_gene_variant', 'downstream_gene_variant', 'TF_binding_site_variant',
                  'regulatory_region_variant', 'miRNA', 'custom', 'sequence_feature',
                  'conserved_intron_variant', 'intron_variant', 'intragenic_variant',
                  'conserved_intergenic_variant', 'intergenic_region', 'coding_sequence_variant',
                  'non_coding_exon_variant', 'nc_transcript_variant', 'gene_variant', 'chromosome']

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

if args.input == '-':
    inhandle = sys.stdin
else:
    inhandle = open(args.input, 'rb')

if args.output is None or args.output == '-':
    outhandle = sys.stdout
else:
    outhandle = open(args.output, 'wb')

with inhandle as ih:
    with outhandle as out:
        for line in ih:
            fields = line.strip().split('\t')
            annot = [j for f in fields[2:] for j in f.split('&')]
            worst = None
            worst_index = len(variant_order)
            for v in annot:
                try:
                    i = variant_order.index(v)
                    if i < worst_index:
                        worst_index = i
                        worst = v
                except ValueError:
                    continue
            if worst is None:
                continue
            worst = re.sub('_variant', '', worst)
            if worst == 'downstream_gene' or worst == 'upstream_gene':
                worst = 'intergenic'
            p = int(fields[1])
            out.write(fields[0] + '\t' + str(p - 1) + '\t' + str(p) +
                      '\t' + worst + '\n')
