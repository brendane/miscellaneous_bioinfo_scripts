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

            rc = False
            if qry_pos[1] < qry_pos[0]:
                rc = True
                tmp_ = qry_pos[0]
                qry_pos[0] = qry_pos[1]
                qry_pos[1] = tmp_

            # Get sequence flanking sequences
            # Assume that reference sequences are each the full
            # sequence of interest
            ref_seq = ref_idx[ref_tag]
            qry_pos_fl = copy.copy(qry_pos)
            qry_pos_fl[0] -= 1000
            qry_pos_fl[1] += 1000
            if qry_pos_fl[0] < 1:
                qry_pos_fl[0] = 1
            if qry_pos_fl[1] > len(qry_idx[qry_tag]):
                qry_pos_fl[1] = len(qry_idx[qry_tag])
            qry_seq = qry_idx[qry_tag][(qry_pos_fl[0]-1):qry_pos_fl[1]]
            if rc:
                qry_fl_1 = qry_idx[qry_tag][(qry_pos_fl[0]-1):(qry_pos[0]-1)]
                qry_fl_0 = qry_idx[qry_tag][qry_pos[1]:(qry_pos_fl[1])]
            else:
                qry_fl_0 = qry_idx[qry_tag][(qry_pos_fl[0]-1):(qry_pos[0]-1)]
                qry_fl_1 = qry_idx[qry_tag][qry_pos[1]:(qry_pos_fl[1])]

            if rc:
                qry_seq.seq = Seq.reverse_complement(qry_seq.seq)
                qry_fl_0.seq = Seq.reverse_complement(qry_fl_0.seq)
                qry_fl_1.seq = Seq.reverse_complement(qry_fl_1.seq)

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

            qry_seql = [qry_seq]
            qry_fl_0l = [qry_fl_0]
            qry_fl_1l = [qry_fl_1]
            if qry_file in ref_tag:
                # TODO
                # This isn't the right way to handle the problem. Instead,
                # take first and last alignment position from this row and
                # the previous rows, then re-extract the sequence.
                #
                # Even that can have issues if there are inversions, but this
                # program is designed for straightforward cases.
                sys.stderr.write('Warning: %s has multiple alignments in %s\n' % 
                                 (ref_tag, qry_file))
                qry_seql = info[ref_tag][qry_file]['qry_seq'] + qry_seql
                qry_fl_0l = info[ref_tag][qry_file]['qry_flank_0'] + qry_fl_0l
                qry_fl_1l = info[ref_tag][qry_file]['qry_flank_1'] + qry_fl_1l
            info[ref_tag][qry_file] = {'query_tag':qry_tag,
                                       'query_pos':qry_pos,
                                       'query_pos_fl':qry_pos_fl,
                                       'qry_seq':qry_seql,
                                       'qry_flank_0':qry_fl_0l,
                                       'qry_flank_1':qry_fl_1l}

# Run alignments on each flanking sequence
for tag in info:
    # It seems like running the alignment on the whole sequence, not
    # just flanking regions, works slightly better.
    #
    # I think this is an issue mainly for genes for which there are
    # multiple nucmer matches, and it is due to a bug.
    """
    for flank in ['0', '1']:
        fd, tmp = tempfile.mkstemp(dir=tempdir)
        with os.fdopen(fd, 'wb') as oh:
            #oh.write('>ref\n')
            #oh.write(str(info[tag]['ref_seq'].seq) + '\n')
            for qf in info[tag]:
                if qf == 'ref_seq':
                    continue
                oh.write('>' + qf + '_flank_' + flank + '\n')
                s = ''.join(map(lambda x: str(x.seq), info[tag][qf]['qry_flank_' + flank]))
                oh.write(s + '\n')

        p = subprocess.Popen(['muscle', '-in', tmp,
                              '-out', tmp + '.aln', '-quiet'])
        retcode = p.wait()
        if retcode != 0:
            raise Exception('muscle did not work')

        # Make consensus sequence and check that region is reasonably
        # similar
        aln = AlignIO.read(tmp + '.aln', 'fasta')

        # Report flanking sequences in a separate file

        # Clean up
        os.unlink(tmp)
        shutil.copyfile(tmp + '.aln',
                    args.output + '/' + tag + '_flank_' + flank + '.fasta')
        os.unlink(tmp + '.aln')
    """

    fd, tmp = tempfile.mkstemp(dir=tempdir)
    with os.fdopen(fd, 'wb') as oh:
        oh.write('>ref\n')
        oh.write(str(info[tag]['ref_seq'].seq) + '\n')
        for qf in info[tag]:
            if qf == 'ref_seq':
                continue
            oh.write('>' + qf + '\n')
            s = ''.join(map(lambda x: str(x.seq), info[tag][qf]['qry_seq']))
            oh.write(s + '\n')

    p = subprocess.Popen(['muscle', '-in', tmp,
                          '-out', tmp + '.aln', '-quiet'])
    retcode = p.wait()
    if retcode != 0:
        raise Exception('muscle did not work')

    # Make consensus sequence and check that region is reasonably
    # similar
    aln = AlignIO.read(tmp + '.aln', 'fasta')
    ## TODO: get start and end of each flanking sequence in the
    ## alignment and make a consensus

    # Clean up
    os.unlink(tmp)
    shutil.copyfile(tmp + '.aln',
                args.output + '/' + tag + '.fasta')
    os.unlink(tmp + '.aln')



## TODO
##  Consensus sequence extraction
##  Manually check results
