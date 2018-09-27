#!/usr/bin/env python2.7
"""
    List variant locations in a bam or cram file from a pileup.
    Best for short sequences.

    Reads mpileup output from stdin and writes to stdout.

    One optional argument: proportion reference allele below which
        the site will be reported.

    Not sure what the incorrect output might be, but the script this
    script is based on had a warning, so I kept it.
"""

import sys

threshold = 1
if len(sys.argv) > 1:
    threshold = float(sys.argv[1])

print '#NOTE: May produce incorrect output under some circumstances'
print 'contig\tbp\tref\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous'
for line in sys.stdin:
    data = line.strip().split('\t')
    contig = data[0]
    bp = data[1]
    bases = data[4].upper()
    ref = data[2].upper()
    
    types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

    i = 0
    reads = 0.
    while i < len(bases):
        base = bases[i].upper()
        if base == '^' or base == '$':
            i += 1
        elif base == '-':
            i += 1
        elif base == '*':
            types['-'] += 1
            reads += 1
        elif base == '+':
            i += 1
            addNum = int(bases[i])
            addSeq = ''
            for a in range(addNum):
                i += 1
                addSeq += bases[i]

            types['+'].append(addSeq)
            reads += 1
        elif base == '.' or base == ',':
            types[ref] += 1
            reads += 1
        else:
            if types.has_key(base):
                types[base] += 1
                reads += 1
            else:
                types['X'].append(base)
                reads += 1

        i += 1

    if types[ref] / reads < threshold:
        adds = '.'
        if len(types['+']) > 0:
            adds = ','.join(types['+'])

        amb = '.'
        if len(types['X']) > 0:
            amb = ','.join(types['X'])


        # For counts
        out = [contig,bp,ref,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb]
        print '\t'.join([str(x) for x in out])
