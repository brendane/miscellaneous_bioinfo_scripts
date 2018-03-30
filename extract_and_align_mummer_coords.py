#!/usr/bin/env python2.7
"""
    Extract and align sequences from MUMmer coords files.

    Note that there are still some issues to be fixed, but they may not
    affect the output if there are simple alignments.
"""

import argparse
import copy
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO, AlignIO, Seq

iupac = {'ACGT':'N',
         'AGT':'D', 'ACG':'V', 'ACT':'H', 'CGT':'B',
         'CT':'Y', 'AG':'R', 'AT':'W', 'CG':'S', 'GT':'K', 'AC':'M'}

def consense(bases):
    if len(bases) == 1 and '-' in bases:
        ret = ''
    elif '-' in bases:
        ret = 'N'
    elif len(bases) == 1:
        ret = list(bases)[0]
    else:
        ret = iupac[''.join(sorted(list(bases)))]
    return ret

def guess_strand(s):
    seq = str(s.seq)
    if seq[:3] == 'ATG' or seq[-3:] in {'TAA', 'TAG', 'TGA'}:
        return '+'
    elif seq[:3] in {'TTA', 'CTA', 'TCA'} or seq[-3:] == 'CAT':
        return '-'
    else:
        return '+'

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--tempdir')
parser.add_argument('--tol', type=float, default=-1.0)
parser.add_argument('--flank', type=int, default=1000)
parser.add_argument('coords', nargs='+')
args = parser.parse_args()

tol = args.tol
flank = args.flank
tempdir = args.tempdir

info = {}
for fname in args.coords:
    with open(fname, 'rb') as ih:
        ref_file, qry_file = ih.readline().strip().split()

        qry_idx = SeqIO.index(qry_file, 'fasta')
        ref_idx = SeqIO.index(ref_file, 'fasta')

        ih.readline(); ih.readline(); ih.readline()
        for line in ih:
            row = line.strip().split('\t')
            ref_tag, qry_tag = row[7:9]
            ref_pos = map(int, row[0:2])
            qry_pos = map(int, row[2:4])
            ref_aln_len, qry_aln_len = map(int, row[4:6])

            # Store information
            # Assume that reference sequences are each the full
            # sequence of interest
            ref_seq = ref_idx[ref_tag]
            if ref_tag not in info:
                info[ref_tag] = {}
                info[ref_tag]['ref_seq'] = ref_seq

            # Check for need to reverse-complement
            rc = False
            strand = '+'
            if qry_pos[1] < qry_pos[0]:
                rc = True
                strand = '-'
                tmp_ = qry_pos[0]
                qry_pos[0] = qry_pos[1]
                qry_pos[1] = tmp_

            # If there is more than one alignment, merge if the
            # strand is the same. If the strand is not the same
            # or the contigs are different, skip with a warning.
            if qry_file in info[ref_tag]:
                if info[ref_tag][qry_file]['strand'] != strand or \
                   info[ref_tag][qry_file]['query_tag'] != qry_tag:
                    sys.stderr.write('Unmergeable alignments for %s from %s\n' %
                                     (ref_tag, qry_file))
                    continue
                qry_pos[0] = min(qry_pos[0], info[ref_tag][qry_file]['query_pos'][0])
                qry_pos[1] = max(qry_pos[1], info[ref_tag][qry_file]['query_pos'][1])

            # Get coordinates of gene + flanking regions. Adjust if
            # initial coordinates are past the boundaries of the contig.
            qry_pos_fl = copy.copy(qry_pos)
            qry_pos_fl[0] -= 1000
            qry_pos_fl[1] += 1000
            if qry_pos_fl[0] < 1:
                qry_pos_fl[0] = 1
            if qry_pos_fl[1] > len(qry_idx[qry_tag]):
                qry_pos_fl[1] = len(qry_idx[qry_tag])

            # Extract the full sequence, including gene and flanking regions
            qry_seq = qry_idx[qry_tag][(qry_pos_fl[0]-1):qry_pos_fl[1]]

            # Extract the just the flanking sequences and do reverse
            # complementing if needed.
            if rc:
                qry_fl_1 = qry_idx[qry_tag][(qry_pos_fl[0]-1):(qry_pos[0]-1)]
                qry_fl_0 = qry_idx[qry_tag][qry_pos[1]:(qry_pos_fl[1])]
                qry_seq.seq = Seq.reverse_complement(qry_seq.seq)
                qry_fl_0.seq = Seq.reverse_complement(qry_fl_0.seq)
                qry_fl_1.seq = Seq.reverse_complement(qry_fl_1.seq)
            else:
                qry_fl_0 = qry_idx[qry_tag][(qry_pos_fl[0]-1):(qry_pos[0]-1)]
                qry_fl_1 = qry_idx[qry_tag][qry_pos[1]:(qry_pos_fl[1])]

            # Do some quality control to make sure the alignments are
            # not weird or missing something
            if tol >= 0:
                if ref_aln_len * (1-tol) > qry_aln_len or \
                   ref_aln_len * (1+tol) < qry_aln_len:
                    sys.stderr.write('Warning: %s may have large indels\n' % ref_tag)
                ref_real_length = len(ref_seq)
                if ref_real_length != ref_aln_len:
                    sys.stderr.write('Warning: %s is not fully aligned\n' % ref_tag)

            # Store all the information
            info[ref_tag][qry_file] = {'query_tag':qry_tag,
                                       'query_pos':qry_pos,
                                       'query_pos_fl':qry_pos_fl,
                                       'qry_seq':qry_seq,
                                       'qry_flank_0':qry_fl_0,
                                       'qry_flank_1':qry_fl_1,
                                       'strand':strand}

# Run alignments on each flanking sequence
with open(args.output + '/summary.tsv', 'wb') as summary:
    summary.write('gene\tn_strains\tstrains\tcds_start\tcds_end\t' +
                  'sequence\tcds_sequence\tflank_1_sequence\t' +
                  'flank_2_sequence\n')
    for tag in info:
        fd, tmp = tempfile.mkstemp(dir=tempdir)
        with os.fdopen(fd, 'wb') as oh:
            oh.write('>ref\n')
            oh.write(str(info[tag]['ref_seq'].seq) + '\n')
            for qf in info[tag]:
                if qf == 'ref_seq':
                    continue
                oh.write('>' + qf + '\n')
                s = str(info[tag][qf]['qry_seq'].seq)
                oh.write(s + '\n')

        p = subprocess.Popen(['muscle', '-in', tmp,
                              '-out', tmp + '.aln', '-quiet'])
        retcode = p.wait()
        if retcode != 0:
            raise Exception('muscle did not work')

        # Make consensus sequence and check that region is reasonably
        # similar.
        # Note that the consensus is the consensus of the target strains,
        # not including the reference.
        aln = AlignIO.read(tmp + '.aln', 'fasta')

        ref_idx = map(lambda x: x.id, aln).index('ref')
        targets_present = [rec.id for rec in aln if rec.id != 'ref']
        consensus = ''
        aln_gene_start = 0
        aln_gene_end = 0
        ref_started = False
        deleted = []
        for bp in xrange(0, len(aln[0])):
            targets = set()
            for rec in aln:
                if rec.id != 'ref':
                    targets.add(rec[bp])
            c =  consense(targets)
            if c == '':
                deleted.append(bp)
            consensus += c
            if not ref_started:
                aln_gene_start = bp
                if aln[ref_idx, bp] != '-':
                    ref_started = True
            elif aln[ref_idx, bp] != '-':
                aln_gene_end = bp

        # Adjust for sites that were gaps in all targets
        before_gene = sum(1 for x in deleted if x <= aln_gene_start)
        in_gene = sum(1 for x in deleted if x <= aln_gene_end)
        consensus_gene_start = aln_gene_start - before_gene
        consensus_gene_end = aln_gene_end - in_gene

        # Make sequences
        consensus_gene = consensus[(consensus_gene_start):(consensus_gene_end+1)]
        consensus_flank_0 = consensus[:consensus_gene_start]
        consensus_flank_1 = consensus[(consensus_gene_end+1):]

        # Write information to output summary file
        output = [tag, len(targets_present),
                  ','.join(sorted(targets_present)),
                  consensus_gene_start+1,
                  consensus_gene_end+1,
                  consensus,
                  consensus_gene,
                  consensus_flank_0,
                  consensus_flank_1 ]
        summary.write('\t'.join(map(str, output)) + '\n')

        # Clean up
        os.unlink(tmp)
        shutil.copyfile(tmp + '.aln',
                    args.output + '/' + tag + '.fasta')
        os.unlink(tmp + '.aln')

# Make summary file for each strain separately
for i, fname in enumerate(args.coords):
    # Get the qry file tag name
    qf = ''
    with open(fname, 'rb') as ih:
        _, qf = ih.readline().strip().split()
    with open(args.output + '/summary_' + str(i) + '.tsv', 'wb') as oh:
        oh.write('# ' + fname + '\n')
        oh.write('gene\tcds\tflank1\tflank2\tfull\n')
        for tag in info:
            if qf not in info[tag]:
                continue
            full_seq = info[tag][qf]['qry_seq']
            f0_seq = info[tag][qf]['qry_flank_0']
            f1_seq = info[tag][qf]['qry_flank_1']
            cds_seq = full_seq[len(f0_seq):(len(full_seq)-len(f1_seq))]
            if guess_strand(cds_seq) == '-':
                full_seq.seq = Seq.reverse_complement(full_seq.seq)
                cds_seq.seq = Seq.reverse_complement(cds_seq.seq)
                tmp = f0_seq
                f0_seq.seq = Seq.reverse_complement(f1_seq.seq)
                f1_seq.seq = Seq.reverse_complement(tmp.seq)
            oh.write(tag + '\t' + str(cds_seq.seq) + '\t' + 
                     str(f0_seq.seq) + '\t' + str(f1_seq.seq) + 
                     '\t' + str(full_seq.seq) + '\n')

