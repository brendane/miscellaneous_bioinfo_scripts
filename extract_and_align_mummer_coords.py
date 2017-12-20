#!/usr/bin/env python2.7
"""
    Extract and align sequences from MUMmer coords files.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO

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

            # Get sequence
            # Assume that reference sequences are each the full
            # sequence of interest
            ref_seq = ref_idx[ref_tag]
            qry_pos_fl = qry_pos
            qry_pos_fl[0] -= 1000
            qry_pos_fl[1] += 1000
            if qry_pos_fl[0] < 1:
                qry_pos_fl[0] = 1
            if qry_pos_fl[1] > len(qry_idx[qry_tag]):
                qry_pos_fl[1] = len(qry_idx[qry_tag])
            qry_seq = qry_idx[qry_tag][(qry_pos_fl[0]-1):qry_pos_fl[1]]

            # Do some quality control to make sure the alignments are
            # not weird or missing something
            if tol >= 0:
                if ref_aln_len * (1-tol) > qry_aln_len or \
                   ref_aln_len * (1+tol) < qry_aln_len:
                    sys.stderr.write('Warning: %s may have large indels\n' % ref_tag)
                ref_real_length = len(ref_seq)
                if ref_real_length != ref_aln_len:
                    sys.stderr.write('Warning: %s is not fully aligned\n' % ref_tag)

            # Store information
            if ref_tag not in info:
                info[ref_tag] = {}
                info[ref_tag]['ref_seq'] = ref_seq
            info[ref_tag][qry_file] = {'query_tag':qry_tag,
                                       'query_pos':qry_pos,
                                       'query_pos_fl':qry_pos_fl,
                                       'qry_seq':qry_seq}

# Run alignments
for tag in info:
    fd, tmp = tempfile.mkstemp(dir=tempdir)
    with os.fdopen(fd, 'wb') as oh:
        oh.write('>ref\n')
        oh.write(str(info[tag]['ref_seq'].seq) + '\n')
        for qf in info[tag]:
            if qf == 'ref_seq':
                continue
            oh.write('>' + qf + '\n')
            oh.write(str(info[tag][qf]['qry_seq'].seq) + '\n')
    p = subprocess.Popen(['muscle', '-in', tmp,
                          '-out', tmp + '.aln', '-quiet'])
    retcode = p.wait()
    if retcode != 0:
        raise Exception('muscle did not work')

    # Make consensus sequence, get alignment start and end coords,
    # and check that flanking regions are reasonably similar

    # Report flanking sequences in a separate file

    # Clean up
    os.unlink(tmp)
    shutil.copyfile(tmp + '.aln',
                args.output + '/' + tag + '.fasta')
    os.unlink(tmp + '.aln')

## TODO
##  Check for sequences that more than one alignment and set aside for
##      manual analysis
##  Consensus sequence extraction
##  Manually check results
