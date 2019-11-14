#!/usr/bin/env python3
"""
    Evaluate a MAF file from mugsy.

    evaluate_mugsy_maf.py --output <output prefix> <input maf file>
        <fai directory> <gff directory> [<comma separated list of strains>]

    Assumes no whitespace in the sequences.
"""

import argparse
import collections
import re
import os
import subprocess
import tempfile

def parse_maf(fname):
    with open(fname, 'rt') as ih:
        record = {'a':{}, 's':{}}
        for line in ih:
            if line.startswith('a'):
                if len(record['s']) > 0:
                    yield record
                record = {'a':{}, 's':{}}
                record['a'] = {x.split('=')[0]:x.split('=')[1] for x in line.strip().split()[1:]}
                record['s'] = {}
            elif line.startswith('s'):
                fields = line.strip().split()[1:]
                record['s'][(fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]))] = fields[5]
            else:
                continue
        yield record


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('maf')
parser.add_argument('fai')
parser.add_argument('gff')
parser.add_argument('strains', nargs='*')
args = parser.parse_args()

n_sites = 0
n_lcbs = 0
n_core_lcbs = 0
lcb_lengths = []
lcb_ungapped_lengths = []
core_lcb_ungapped_lengths = []
lcb_strain_ungapped_lengths = collections.defaultdict(list)
core_lcb_strain_ungapped_lengths = collections.defaultdict(list)
n_multiallelic_sites = 0
lcb_boundaries = collections.defaultdict(list)
lcb_bed = collections.defaultdict(list)
gap_lengths = []
n_multiallelic_sites = 0
all_strains = set()

if args.strains is None:
    n_strains = len(os.listdir(args.fai))
else:
    n_strains = len(args.strains[0].strip().split(','))

for rec in parse_maf(args.maf):
    names = list(rec['s'].keys())
    strains = set(x[0].split('.', 1)[0] for x in names)
    all_strains.update(strains)
    core = False
    if len(strains) == n_strains:
        core = True
        n_core_lcbs += 1
    n_lcbs += 1
    l = len(rec['s'][names[0]])
    lcb_lengths.append(l)
    u = 0
    for i in range(l):
        gap = False
        alleles = set()
        for n in names:
            if rec['s'][n][i] == '-':
                gap = True
            else:
                alleles.add(rec['s'][n][i].upper())
        if not gap: u += 1
        if len(alleles) > 2:
            n_multiallelic_sites += 1
    lcb_ungapped_lengths.append(u)
    if core:
        core_lcb_ungapped_lengths.append(u)
    for n in names:
        gap_lengths += [len(x) for x in re.split('[ATCGNatcgn]',
                                                 rec['s'][n]) if len(x) > 0]
        ll = len(rec['s'][n].replace('-', ''))
        strain = n[0].split('.', 1)[0]
        replicon = n[0].split('.', 1)[1]
        lcb_strain_ungapped_lengths[strain].append(ll)
        if core:
            core_lcb_strain_ungapped_lengths[strain].append(ll)
        if n[3] == '+':
            s = n[1] + n[2] - 1
            e = s + 1
            ss = n[1]
            ee = n[1] + n[2]
        else:
            ## NOTE: Revised 2019-11-14
            s = n[4] - n[1] - n[2] - 1
            e = s + 1
            ee = n[4] - n[1]
            ss = n[4] - n[1] - n[2] - 1
        lcb_boundaries[strain].append((replicon, s, e))
        lcb_bed[strain].append((replicon, ss, ee))

n_sites = sum(lcb_lengths)  # Total length of all blocks
n_ungapped_sites = sum(lcb_ungapped_lengths) # Length of all blocks at ungapped sites
n_core_sites = sum(core_lcb_ungapped_lengths) # Number of bp found in every strain

t = 0
for l in sorted(lcb_lengths, reverse=True):
    t += l
    if t > n_sites / 2:
        n50 = l # N50 for the entire set of LCBs
        break

n_sites_strains = {}  # Number of sites included in the file for each strain
n_core_lcb_sites_strains = {} # Number of sites included in core LCBs for each strain (not necessarily ungapped in all other strains)
n50_strains = {}
for strain in all_strains:
    n_sites_strains[strain] = sum(lcb_strain_ungapped_lengths[strain])
    n_core_lcb_sites_strains[strain] = sum(core_lcb_strain_ungapped_lengths[strain])
    t = 0
    for l in sorted(lcb_strain_ungapped_lengths[strain], reverse=True):
        t += l
        if t > n_sites_strains[strain] / 2:
            n50_strains[strain] = l # N50 for a particular strain
            break

strain_total_length = {}
strain_cut_genes = {}
for strain in all_strains:
    t = 0
    with open(args.fai + '/' + strain + '.fasta.fai', 'rt') as ih:
        for line in ih:
            t += int(line.strip().split()[1])
    strain_total_length[strain] = t
    tmp = tempfile.mkstemp(dir='.')[1]
    with open(tmp, 'wt') as bed:
        for r, s, e in lcb_boundaries[strain]:
            bed.write(r + '\t' + str(s) + '\t' + str(e) + '\n')
    p = subprocess.Popen(['bedtools', 'intersect', '-a',
                          tmp, '-wb', '-b',
                          args.gff + '/' + strain + '.gff3'],
                         stdout=subprocess.PIPE)
    x = p.communicate()[0].decode('utf-8').split('\n')
    os.unlink(tmp)
    #with open(args.output + '.lcb_boundary_features.' + strain + '.bed', 'wt') as oh:
    #    oh.write('\n'.join(x) + '\n')
    strain_cut_genes[strain] = len(x)

with open(args.output + '.stats.tsv', 'wt') as oh:
    oh.write('ALL\tn_sites\t' + str(n_sites) + '\t' +
             'Total length of all LCBs\n')
    oh.write('ALL\tn_lcbs\t' + str(n_lcbs) + '\t' +
             'Number of LCBs (alignment chunks)\n')
    oh.write('ALL\tn_core_lcbs\t' + str(n_core_lcbs) + '\t' +
             'Number of LCBs with all strains present\n')
    oh.write('ALL\tn_multiallelic_sites\t' + str(n_multiallelic_sites) +
             '\tNumber of LCB sites with > 2 non-gap alleles\n')
    oh.write('ALL\tproportion_multiallelic_sites\t' +
             str(n_multiallelic_sites / n_sites) +
             '\tProportion of LCB sites with > 2 non-gap alleles\n')
    oh.write('ALL\tcore_bp\t' + str(n_core_sites) + '\t' +
             'Number of sites with all strains present\n')
    oh.write('ALL\tn_ungapped_sites\t' + str(n_ungapped_sites) +
             '\tNumber of alignment columns with no gaps, including singleton LCBs\n')
    oh.write('ALL\tN50\t' + str(n50) + '\t' +
             'N50 for LCBs\n')
    oh.write('ALL\tproportion_core\t' +
             str(n_core_sites / n_sites) +
             '\tNumber of core bp divided by total\n')
    for strain in all_strains:
        try:
            oh.write(strain + '\tN50\t' + str(n50_strains[strain]) +
                     '\tN50 for %s\n' % strain)
        except:
            import pdb; pdb.set_trace()
        oh.write(strain + '\tn_sites\t' + str(n_sites_strains[strain]) +
                 '\tNumber of sites included from %s\n' % strain)
        oh.write(strain + '\tproportion_included\t' +
                 str(n_sites_strains[strain]/strain_total_length[strain]) +
                 '\tProportion of %s included in the maf file\n' % strain)
        oh.write(strain + '\tn_core_lcb_sites\t' +
                 str(n_core_lcb_sites_strains[strain]) +
                 '\tNumber of sites in %s found in core LCBs\n' % strain)
        oh.write(strain + '\tproportion_core_lcb_sites\t' +
                 str(n_core_lcb_sites_strains[strain] / strain_total_length[strain]) +
                 '\tProportion of sites in %s found in core LCBs\n' % strain)
        oh.write(strain + '\tproportion_core_sites\t' +
                 str(n_core_sites / strain_total_length[strain]) +
                 '\tProportion of sites in %s found in core bp\n' % strain)
        oh.write(strain + '\tn_cut_genes\t' +
                 str(strain_cut_genes[strain]) + '\t' +
                 'Number of LCB boundaries that are in the middle of genes in %s\n' % strain)
    for g in sorted(gap_lengths):
        oh.write('GAP\t\t' + str(g) + '\t\n')
