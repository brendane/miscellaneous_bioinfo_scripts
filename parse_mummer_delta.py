#!/usr/bin/env python3
## TODO:
## Not done -- needs revision and testing

import sys

def process_indels(indels, r, q, rs, re, qs, qe):
    rp = rs - 1
    qp = qs - 1
    qdir = 1
    if qs > qe:
        qdir = -1
        qp = qs
    rdir = 1
    if rs > re:
        rdir = -1
        rp = rs
    blocks = []
    for i in indels:
        old_rp = rp
        old_qp = qp
        if i < 0:
            rp += (-i - 1) * rdir
            qp += -i * qdir
        elif i > 0:
            rp += i * rdir
            qp += (i - 1) * qdir
        else:
            rp = re
            if rdir == -1: rp = re -1
            qp = qe
            if qdir == -1: qp = qe - 1
        blocks.append((r, old_rp, rp, q, old_qp, qp))
    import pdb; pdb.set_trace()
    return blocks

alignments = []
with open(sys.argv[1], 'rt') as delta_handle:
    seq_files = delta_handle.readline().strip().split()
    data_type = delta_handle.readline().strip()
    newblock = False
    indels = []
    while True:
        try:
            line = delta_handle.readline().strip()
            if line == '':
                break
        except EOFError:
            break
        if line.startswith('>'):
            if len(indels):
                alignments.append(process_indels(indels, ref, qry, refstart,
                                                 refend, qrystart, qryend))
            ref, qry = line.split()[:2]; ref = ref[1:]
            reflen, qrylen = (int(x) for x in line.split()[2:])
            newblock = True
            indels = []
        elif newblock:
            refstart, refend, qrystart, qryend, err, simerr, stops = (int(x) for x in line.split())
            newblock = False
        else:
            indels.append(int(line))
    alignments.append(process_indels(indels, ref, qry, refstart,
                                     refend, qrystart, qryend))

