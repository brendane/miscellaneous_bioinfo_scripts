#!/usr/bin/env python2.7
"""
    Simple script from http://blog.nextgenetics.net/?e=56

    Often seems to give slightly different counts than I get by hand
    from the mpileup output - probably not able to handle certain cases.

"""

import sys

inFile = sys.stdin

print 'bp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous'
for line in inFile:
        data = line.strip().split('\t')
        bp = data[1]
        bases = data[4].upper()
        ref = data[2].upper()
        
        types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

        i = 0
        while i < len(bases):
                base = bases[i]
                if base == '^' or base == '$':
                        i += 1
                elif base == '-':
                        i += 1
                elif base == '*':
                        types['-'] += 1
                elif base == '+':
                        i += 1
                        addNum = int(bases[i])
                        addSeq = ''
                        for a in range(addNum):
                                i += 1
                                addSeq += bases[i]

                        types['+'].append(addSeq)
                elif base == '.' or base == ',':
                        types[ref] += 1
                else:
                        if types.has_key(base):
                                types[base] += 1
                        else:
                                types['X'].append(base)

                i += 1

        adds = '.'
        if len(types['+']) > 0:
                adds = ','.join(types['+'])

        amb = '.'
        if len(types['X']) > 0:
                amb = ','.join(types['X'])


        print 'NOTE: May produce incorrect output under some circumstances'
        # For counts
        out = [bp,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb]
        print '\t'.join([str(x) for x in out])

        # For frequencies
        total = float(types['A']+types['G']+types['C']+types['T']+
                      types['-']+len(types['+']))
        out = [bp,
               round(types['A']/total, 3),
               round(types['G']/total, 3),
               round(types['C']/total, 3),
               round(types['T']/total, 3),
               round(types['-']/total, 3),
               round(len(types['+'])/total, 3),
               adds,amb]
        print '\t'.join([str(x) for x in out])
