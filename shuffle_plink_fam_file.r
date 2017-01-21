#!/usr/bin/env Rscript

cargs = commandArgs(trailingOnly=TRUE)

infile = cargs[1]
outfile = cargs[2]

x = read.table(infile)
x[, 6] = x[sample(1:nrow(x), nrow(x), FALSE), 6]

write.table(x, outfile, sep=' ', col.names=FALSE, row.names=FALSE,
            quote=FALSE)
