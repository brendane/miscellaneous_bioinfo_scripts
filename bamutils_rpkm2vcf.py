#!/usr/bin/env python
"""
    Convert bamutils RPKM output to VCF.

    bamutils_rpkm2vcf.py --greater <log2 threshold for increase>
        --ambig <log2 threshold for ambiguous>
        --absent <log2 threshold for absent>
        --pad <padding for position>
        <input directory> <reference strain name>

    Output is to stdout.

    The input directory should have files named <strain>.<replicon>.tsv.
    
    Assumes input is in 0-based coordinates.

    --absent < --ambig < --greater
"""

#==============================================================================#

import argparse
import csv
import math
import os
import os.path as osp
import sys

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--absent', type=float)
parser.add_argument('--ambig', type=float)
parser.add_argument('--greater', type=float)
parser.add_argument('--pad', type=int, default=0)
parser.add_argument('indir')
parser.add_argument('refstrain')
args = parser.parse_args()

absent_threshold = args.absent
ambig_threshold = args.ambig
greater_threshold = args.greater
pad = args.pad
ref_strain = args.refstrain
indir = args.indir

if absent_threshold >= ambig_threshold:
    raise Exception('--absent must be < --ambig')
if ambig_threshold >= greater_threshold:
    raise Exception('--greater must be > --ambig')

strains_ = set()
replicons = set()
nfiles = 0
for fname in os.listdir(indir):
    if fname.endswith('.tsv'):
        strain, replicon, _ = fname.split('.')
        strains_.add(strain)
        replicons.add(replicon)
        nfiles += 1
if nfiles != len(strains_) * len(replicons):
    raise Exception('Some strain/replicons missing or extra files')

strains = sorted(strains_)

ref_data = {}
for replicon in replicons:
    fname = osp.join(indir, ref_strain + '.' + replicon + '.tsv')
    with open(fname) as handle:
        lp = handle.tell()
        while handle.readline().startswith('#'):
            lp = handle.tell()
            continue
        handle.seek(lp)
        rdr = csv.DictReader(handle, delimiter='\t')
        for row in rdr:
            ref_data[(row['name'])] = (float(row['RPKM']), row['chrom'],
                                       int(row['start']), int(row['end']))

strain_data = {}
for strain in strains:
    strain_data[strain] = {}
    for replicon in replicons:
        fname = osp.join(indir, strain + '.' + replicon + '.tsv')
        with open(fname) as handle:
            lp = handle.tell()
            while handle.readline().startswith('#'):
                lp = handle.tell()
                continue
            handle.seek(lp)
            rdr = csv.DictReader(handle, delimiter='\t')
            for row in rdr:
                rpkm = float(row['RPKM'])
                gene = row['name']
                if rpkm == 0.:
                    ratio = absent_threshold - 1.
                else:
                    try:
                        ratio = math.log(rpkm / ref_data[gene][0], 2)
                    except ZeroDivisionError:
                        # Hopefully, this never happens, but just in case...
                        ratio = float('nan')
                if ratio > greater_threshold:
                    call = 'G'
                elif ambig_threshold < ratio <= greater_threshold:
                    call = 'T'
                elif ratio < absent_threshold:
                    call = 'A'
                else:
                    call = 'N'
                strain_data[strain][gene] = call


sys.stdout.write('##fileformat=VCFv4.1\n')
sys.stdout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
sys.stdout.write('\t' + '\t'.join(strains) + '\n')

for gene in ref_data:
    chrom = ref_data[gene][1]
    pos = ref_data[gene][2] + pad + 1
    var_id = gene
    gts = []
    for strain in strains:
        try:
            gts.append(strain_data[strain][gene])
        except:
            import pdb; pdb.set_trace()
    unique_gts = set(gts) - set('N')
    #if len(unique_gts) < 2:
    #    continue

    ref_gt = 'T'
    alt_gts = list(unique_gts - set(['T', 'N']))
    if len(alt_gts) == 0:
        alt_gts = ['.']

    sys.stdout.write(chrom + '\t' + str(pos) + '\t' + var_id + '\tT\t')
    sys.stdout.write(','.join(alt_gts) + '\t.\t.\t.\tGT')
    for gt in gts:
        if gt == 'N':
            sys.stdout.write('\t./.')
        elif gt == ref_gt:
            sys.stdout.write('\t0/0')
        else:
            g = str(alt_gts.index(gt) + 1)
            sys.stdout.write('\t' + g + '/' + g)

    sys.stdout.write('\n')
