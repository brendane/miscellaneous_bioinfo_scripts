#!/usr/bin/env Rscript
#
# Given a tsv file with strain frequencies, choose the top
# and bottom strains and output to two files.
#
# The input file should have a column labeled "pool" and one
# additional column for every strain.
#
# choose_top_bottom_strains.r <N> <pool> <infile> <out prefix>
#
# Note that in the case of ties, more than N strains can be
# included.
#

cargs = commandArgs(trailingOnly=TRUE)
n = as.numeric(cargs[1])
pool = cargs[2]
infile = cargs[3]
outpre = cargs[4]

x = read.csv(infile, sep='\t', header=TRUE, as.is=TRUE,
             check.names=FALSE)
y = as.numeric(x[x[, 'pool'] == pool, -1])
names(y) = colnames(x)[-1]

ord = order(y)
low = y[ord[n]]
high = y[ord[length(y)-n+1]]

cat(paste(names(y)[y >= high], collapse='\n'), '\n', sep='',
    file=paste0(outpre, '.top.txt'))
cat(paste(names(y)[y <= low], collapse='\n'), '\n', sep='',
    file=paste0(outpre, '.bottom.txt'))

cat(low, '\n')
cat(high, '\n')
