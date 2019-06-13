#!/usr/bin/env python3
"""
    Using the tsv output from parse_mummer_delta.py,
    combine orthologous sites (only 1-1 sites should be in the 
    input).

    align_vcfs_from_genome_alignment.py --output <output vcf> 
        <alignment file> <reference genome 2> <vcf1> <vcf2>

    Note that this can only handle simple SNPs. Remove indels and MNPs
    from both VCF files before running. Assumes haploid genomes.

    vcf1 should be called against the genome given as the reference to
    blat (listed first in command), which will be in the 3rd and 4th
    columns. vcf1 should also be sorted by position and chromosome.

    Requires the pysam and Biopython modules.
"""

import argparse
import csv

from Bio import SeqIO, Seq
import pysam

def strand_match(b, strand):
    if strand == '+':
        return b
    else:
        return Seq.complement(b)

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('aln')
parser.add_argument('ref2')
parser.add_argument('vcf1')
parser.add_argument('vcf2')
args = parser.parse_args()

## Read alignment file
aln = {}
with open(args.aln, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        aln[(row[2], int(row[3]))] = (row[0], int(row[1]), row[4])

## Read reference genome
ref2 = {}
for rec in SeqIO.parse(args.ref2, 'fasta'):
    ref2[rec.id] = str(rec.seq)

## Get vcf files
vf2 = pysam.VariantFile(args.vcf2)
vf1 = pysam.VariantFile(args.vcf1)
samples1 = [s for s in vf1.header.samples]
samples2 = [s for s in vf2.header.samples]

## Make output handle
vf_out = pysam.VariantFile(args.output, 'w', header=vf1.header)
for s in samples2:
    vf_out.header.add_sample(s)
vf_out.close()

## Loop over VCF1 and write matching records
old_pos1 = ('.', 0)
with open(args.output, 'a') as oh:
    for rec in vf1:
        ## Skip MNPs and indels
        if max(map(len, rec.alleles)) > 1:
            continue
        pos1 = (rec.chrom, rec.pos)

        ## Get matching position; if no matching position, skip this position
        try:
            c2, p2, strand = aln[pos1] 
        except KeyError:
            continue
        pos2 = (c2, p2)

        ## TODO -- Check for positions between this one and the last
        ## position in the VCF file. If there are fixed differences
        ## between the references, these should be added.
        ## Problem: We don't know whether these are covered or genotyped
        ## in either species if they aren't in the VCF. Need a VCF with
        ## all positions, or a different approach.
        #if old_pos1[0] != '.' and old_pos1[0] != pos1[0]:
        #    ## Check to end of previous chromosome
        #    pass
        #elif old_pos1[1] != pos1[1]:
        #    ## From previous position to here
        #    for p in range(old_pos1[1]+1, pos1[1]):
        #        pos1_ = (rec.chrom, p)
        #        try:
        #            cc2, pp2, st = aln[pos1] 
        #        except KeyError:
        #            continue
        #        pos2_ = (cc2, pp2)

        ## Another problem: This script does not deal with variant sites
        ## in vcf2 unless they have a match in vcf1.


        ## Get alleles and genotypes
        genotypes = []
        for s in samples1:
            gt = rec.samples[s].alleles[0]
            if gt is None:
                gt = 'N'
            genotypes.append(gt)
        vcf2_match = [r for r in vf2.fetch(pos2[0], pos2[1]-1, pos2[1])]
        if len(vcf2_match) == 0:
            ## Get reference allele from reference genome
            gt = strand_match(ref2[pos2[0]][pos2[1]-1], strand)
            for s in samples2:
                genotypes.append(gt)
        else:
            for s in samples2:
                gt = vcf2_match[0].samples[s].alleles[0]
                if gt is None:
                    gt = 'N'
                genotypes.append(strand_match(gt, strand))
            pass

        ## Construct a new record and write output; handle this with
        ## text to keep it simple.
        alleles = set(genotypes) - set('N')
        if len(alleles) == 1:
            continue
        alts = list(alleles - set(rec.ref))
        allele_idx = [rec.ref] + alts
        out = [pos1[0], str(pos1[1]), rec.id, rec.ref, ','.join(alts),
               '.', '.', '.', 'GT',
               '\t'.join(str(allele_idx.index(gt)) if gt != 'N' else '.' for gt in genotypes)]
        oh.write('\t'.join(out) + '\n')
