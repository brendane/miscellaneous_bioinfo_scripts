#!/usr/bin/env python3
"""
    Extract a region around a reference genome position given the
    output from show-aligns in MUMmer and print the region as two
    fasta records.

    extract_ref_region_from_show_aligns.py --output <output prefix>
        <delta file> <variant bed file> <alignment bed file>

    Requires that show-aligns is available on the PATH and mafft.

    The bed file should have positions of interest in the reference
    genome, and the fourth column should be query genome coordinates
    to run show-aligns on: contig:start-end.
"""

import argparse
import io
import os
import subprocess
import sys
import tempfile

from Bio import AlignIO

def process_aln_txt(a):
    ret = []
    r_begin = None
    r_end = None
    q_begin = None
    q_end = None
    in_aln = False
    r_seq = ""
    q_seq = ""
    for line in a:
        line = line.strip()
        if line.startswith('-- BEGIN alignment'):
            in_aln = True
            fields = line.split()
            r_begin = int(fields[5])
            r_end = int(fields[7])
            q_begin = int(fields[10])
            q_end = int(fields[12])
        elif line.startswith('--   END alignment'):
            in_aln = False
            ret.append(((r_begin, r_end, q_begin, q_end), r_seq, q_seq))
            r_seq = ''
            q_seq = ''
        elif in_aln:
            if line == '':
                continue
            if '|' in line or '^' in line:
                continue
            s = ''.join(line.split()[1:]).upper()
            if len(r_seq) > len(q_seq):
                q_seq += s
            else:
                r_seq += s
        else:
            continue
    return ret

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('delta')
parser.add_argument('vars')
parser.add_argument('matches')
args = parser.parse_args()

    
delta_file = args.delta
var_file = args.vars
bed_file = args.matches

variants = []
with open(var_file, 'rt') as ih:
    for line in ih:
        variants.append(tuple(line.strip().split()[:3]))

bed_entries = []
with open(bed_file, 'rt') as ih:
    for line in ih:
        bed_entries.append(tuple(line.strip().split()[:4]))


for entry in bed_entries:
    ref_contig = entry[0]
    r_start, r_end = map(int, entry[1:3])
    query_contig = entry[3].split(':')[0]
    q_start, q_end = map(int, entry[3].split(':')[1].split('-'))
    for v in variants:
        if int(v[1]) >= r_start and int(v[2]) <= r_end:
            var_contig = v[0]
            var_start = int(v[1])
            var_end = int(v[2])

    s_a_process = subprocess.Popen(['show-aligns', '-w', '1000000000', delta_file, ref_contig, query_contig],
                                   stdout=subprocess.PIPE)
    aln_txt = s_a_process.communicate()[0].decode('utf-8')
    alns = process_aln_txt(aln_txt.split('\n'))
    for (s, e, _, _), r_seq, q_seq in alns:
        ## If this alignment contains the target
        if s > r_start and e <= r_end:
            k = 0
            ## Find the target position in the alignment
            s_pos, s_pos_ref = None, None
            e_pos, e_pos_ref  = None, None
            for i, b in enumerate(r_seq):
                if b != '-' and b != '.':
                    k += 1
                if k + r_start - 1 == var_start:
                    s_pos = i
                    s_pos_ref = k - 1
                if k + r_start == var_end:
                    e_pos = i + 1
                    e_pos_ref = k
                    break
            # Re-align 300 bp on either side using MAFFT
            if s_pos < 300:
                slice_start = 0
                s_pos_in_slice = s_pos_ref
                e_pos_in_slice = e_pos_ref
            else:
                slice_start = s_pos - 300
                d = e_pos - s_pos; dd = e_pos_ref - s_pos_ref
                s_pos_in_slice = 300 - r_seq[slice_start:(slice_start+300)].count('.')
                e_pos_in_slice = 300 - r_seq[slice_start:(slice_start+301+d)].count('.') + dd
            if len(r_seq) - e_pos < 300:
                slice_end = len(r_seq)
            else:
                slice_end = e_pos + 300
            r_seq = r_seq[slice_start:slice_end].replace('.', '').replace('-', '')
            q_seq = q_seq[slice_start:slice_end].replace('.', '').replace('-', '')
            fd, alnfile = tempfile.mkstemp(dir='.')
            os.close(fd)
            with open(alnfile, 'wt') as handle:
                handle.write('>ref\n' + r_seq + '\n' + '>qry\n' + q_seq)
            fasta_alignment = subprocess.Popen(['linsi', alnfile],
                                               stdout=subprocess.PIPE)
            re_alignment = AlignIO.read(io.StringIO(fasta_alignment.communicate()[0].decode('utf-8')),
                                        'fasta')
            k = 0
            r_allele = ''; q_allele = ''
            begin = 0
            end = 0
            for i, b in enumerate(re_alignment[0]):
                if k >= s_pos_in_slice and k < e_pos_in_slice:
                    if begin == 0: begin = i
                    if end == 0: end = i+1
                    r_allele += b
                    q_allele += re_alignment[1][i]
                if b != '-':
                    k += 1
            beg = begin - 15
            ed = end + 15
            if beg < 0: beg = 0
            if ed > len(re_alignment[0]): ed = len(re_alignment[0])
            print('\t'.join(map(str, [var_contig, var_start, var_end, r_allele, q_allele,
                                      re_alignment[0][beg:ed].seq,
                                      re_alignment[1][beg:ed].seq])))
            os.unlink(alnfile)
