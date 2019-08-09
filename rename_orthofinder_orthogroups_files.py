#!/usr/bin/env python3
"""
    Change gene names in the Orthogroups output files from OrthoFinder
    v2.2.

    rename_orthofinder_orthogroups_files.py --output-prefix <output prefix>
        <new SequenceIDs.txt> <SpeciesIDs.txt> <MCL file>

    This script creates new files and only makes Orthgroups.csv, Orthogroups.txt,
    and Orthogroups_UnassignedGenes.csv.
"""

import argparse
import collections

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output-prefix')
parser.add_argument('si')
parser.add_argument('spp')
parser.add_argument('mcl')
args = parser.parse_args()


## Use the correct SequenceIDs.txt file to get the names of genes
ids = {}
with open(args.si, 'rt') as ih:
    for line in ih:
        ID, lt = line.split(': ')
        ids[ID] = lt


## Use the SpeciesIDs.txt file to get the names of genes
spp_ids = {}
with open(args.spp, 'rt') as ih:
    for line in ih:
        ID, sp = line.split(': ')
        spp_ids[ID] = sp


## Read in the MCL file to get the assignment of genes to orthogroups
orthogroups = collections.defaultdict(list)
done = set()
begin = False
og = 0
with open(args.mcl, 'rt') as ih:
    for line in ih:
        if begin:
            fields = [f for f in line.split() if f != '$']
            if line.startswith(')'):
                continue
            elif not line.startswith(' '):
                og = 'OG' + fields[0].zfill(7)
                orthogroups[og].extend(fields[1:])
                if len(done.intersection(fields[1:])) > 0:
                    raise Exception('Found repeated genes in MCL file: %s' % ' '.join(fields))
                done.update(fields[1:])
            else:
                orthogroups[og].extend(fields)
                if len(done.intersection(fields)) > 0:
                    raise Exception('Found repeated genes in MCL file: %s' % ' '.join(fields))
                done.update(fields)
        elif begin and line.startswith('begin'):
            begin = True


## Make new Orthgroups.csv file and Orthgroups_UnassignedGenes.csv file
with open(args.output_prefix + '.csv', 'wt') as oh:
    with open(args.output_prefix + '_UnassignedGenes.csv', 'wt') as ohu:
        for i in range(len(spp_ids)):
            oh.write('\t' + spp_ids[str(i)])
            ohu.write('\t' + spp_ids[str(i)])
        oh.write('\n')
        ohu.write('\n')
        for og in sorted(orthogroups):
            genes = orthogroups[og]
            if len(orthogroups[og]) > 1:
                oh.write(og)
            else:
                ohu.write(og)
            for i in range(len(spp_ids)):
                g = ', '.join(gg for gg in genes if gg.startswith(str(i) + '_'))
                if len(orthogroups[og]) > 1:
                    oh.write('\t' + g)
                else:
                    ohu.write('\t' + g)
            if len(orthogroups[og]) > 1:
                oh.write('\n')
            else:
                ohu.write('\n')


## Make new Orthogroups.txt file
with open(args.output_prefix + '.txt', 'wt') as oh:
    for og in sorted(orthogroups):
        genes = orthogroups[og]
        oh.write(og + ': ' + ' '.join(sorted(orthogroups[og])) + '\n')
