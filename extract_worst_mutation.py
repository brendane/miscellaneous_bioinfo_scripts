#!/usr/bin/env python2
"""
    From a tab-delimited file from snpEff, extract the most deleterious
    mutation type for each variant.

    Uses the mutation order from:
    http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf

    Changes upstream_gene_variant and downstream_gene_variant to
    intergenic and removes "_variant" from all the names.

    extract_worst_mutation.py --output <output file> <input file>

        --pos-grouping: Group variants by position; if multiple lines
            have the position, collapse them into one variant.
        --rs: third field in input is the rs value; add to last column
            of output

    If input or output is given as "-", read from / write to stdin/stdout.

    The input file should have the annotations in the third column.
"""

#==============================================================================#

import argparse
import csv
import itertools
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

def pass_through(it):
    for i in it:
        yield [(i[0], i[1]), [i]] # List form to match group_by

def chr_pos_grouper(it):
    return itertools.groupby(it, lambda x: (x[0], x[1]))

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--pos-grouping', default=False, action='store_true')
parser.add_argument('--rs', default=False, action='store_true')
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

if args.pos_grouping:
    iterator = chr_pos_grouper
else:
    iterator = pass_through

with inhandle as ih:
    with outhandle as out:
        rdr = csv.reader(ih, delimiter='\t')
        for (chrom, p), rows in iterator(rdr):
            rows = [r for r in rows]
            if args.rs:
                rs = rows[0][2]
                annot_idx = 3
            else:
                annot_idx = 2
            annot = [j for row in rows for f in row[annot_idx:] for j in f.split('&')]
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
            p = int(p)
            if args.rs:
                out.write(chrom + '\t' + str(p - 1) + '\t' + str(p) +
                          '\t' + worst + '\t' + rs + '\n')
            else:
                out.write(chrom + '\t' + str(p - 1) + '\t' + str(p) +
                          '\t' + worst + '\n')
