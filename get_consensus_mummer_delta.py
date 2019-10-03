#!/usr/bin/env python3
"""
    Extract a consensus sequence from a MUMmer delta file (for nucleotides).

    get_consensus_mummer_delta.py <delta file>

    Requires that show-aligns is available on the PATH.
"""

import subprocess
import sys

iupac = {'ACGTN':'N',
         'AGT':'D', 'ACG':'V', 'ACT':'H', 'CGT':'B',
         'CT':'Y', 'AG':'R', 'AT':'W', 'CG':'S', 'GT':'K', 'AC':'M',
         'A':'A', 'C':'C', 'G':'G', 'T':'T'}

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
    
def consensus(s1, s2):
    ret = ''
    for b1, b2 in zip(s1, s2):
        if b1 == b2:
            ret += b1
        else:
            try:
                ret += iupac[(min(b1, b2) + max(b1, b2))]
            except KeyError:
                if '.' in b1:
                    ret += b2
                elif '.' in b2:
                    ret += b1
                else:
                    raise Exception('Unknown nucleotide %s %s' % (b1, b2))
    return ret

delta_file = sys.argv[1]

contigs = set()
with open(sys.argv[1], 'rt') as ih:
    for line in ih:
        if line.startswith('>'):
            r, q, _, _ = line.strip().split()
            contigs.add((r[1:], q))

for r, q in contigs:
    s_a_process = subprocess.Popen(['show-aligns', delta_file, r, q],
                                   stdout=subprocess.PIPE)
    aln_txt = s_a_process.communicate()[0].decode('utf-8')
    alns = process_aln_txt(aln_txt.split('\n'))
    for a in alns:
        sys.stdout.write('>' + r + '_' + str(a[0][0]) + '-' + str(a[0][1]) +
                         '__' + q + '_' + str(a[0][2]) + '-' + str(a[0][3]) + '\n')
        sys.stdout.write(consensus(a[1], a[2]) + '\n')
