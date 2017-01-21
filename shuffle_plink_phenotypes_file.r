#!/usr/bin/env Rscript

cargs = commandArgs(trailingOnly=TRUE)

infile = cargs[1]
outfile = cargs[2]

x = read.table(infile)
x[, 3] = x[sample(1:nrow(x), nrow(x), FALSE), 3]

write.table(x, outfile, sep='\t', col.names=FALSE, row.names=FALSE,
            quote=FALSE)
