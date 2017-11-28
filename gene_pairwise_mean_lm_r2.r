#!/usr/bin/env Rscript
#
# Calculate mean pairwise R^2 between variants in genes using a linear
# regression of all the variants in one gene on each variant in the
# other gene.
#
# Will spit out a bunch of warnings about "perfect fit" -- ignore these.
#

library(data.table)
library(dplyr)
library(dtplyr)

maf = function(g) {
    r = sum(g, na.rm=TRUE) / sum(!is.na(g))
    min(r, 1-r)
}

argv = commandArgs(trailingOnly=TRUE)
bedfile = argv[1]
genofile = argv[2]
min_maf = as.numeric(argv[3])
chunk = as.numeric(argv[4])
nchunks = as.numeric(argv[5])

genos = fread(genofile, header=TRUE, showProgress=FALSE)
mafs = apply(as.matrix(genos[, -(1:3)]), 1, maf)
genos = genos[mafs >= min_maf, ]

g = t(as.matrix(genos[, -(1:3), with=FALSE]))

genes = fread(bedfile, header=FALSE)
gene_ids = gsub('ID=(.+?);.+', '\\1', genes[[9]])
chunk_size = ceiling(nrow(genes) / nchunks)
first = (chunk - 1) * chunk_size + 1
last = chunk * chunk_size
if(last > nrow(genes)) { last = nrow(genes) }
## Note that issues with rounding to the nearest whole number mean
## that sometimes the last few chunks will be empty and cause an
## error.

result = matrix(nrow=(last-first+1), ncol=nrow(genes))
colnames(result) = gene_ids
rownames(result) = gene_ids[first:last]

k = 0
for(i in first:last) {
    for(j in 1:nrow(genes)) {
        if(!is.na(result[i-first+1, j])) {
            next
        }
        g1_vars = genos[['contig']] == genes[i][[1]] &
            genos[['pos']] >= genes[i][[4]] &
            genos[['pos']] <= genes[i][[5]]
        g2_vars = genos[['contig']] == genes[j][[1]] &
            genos[['pos']] >= genes[j][[4]] &
            genos[['pos']] <= genes[j][[5]]
        if(sum(g1_vars) && sum(g2_vars)) {
            g1 = g[, g1_vars, drop=FALSE]
            g2 = g[, g2_vars, drop=FALSE]
            r2_1 = 0
            r2_n1 = 0
            r2_2 = 0
            r2_n2 = 0
            for(l in 1:ncol(g2)) {
                r2_ = summary(lm(g2[, l] ~ g1))$adj.r.squared
                if(!is.na(r2_)) {
                    r2_1 = r2_1 + r2_
                    r2_n1 = r2_n1 + 1
                }
            }
            for(l in 1:ncol(g1)) {
                r2_ = summary(lm(g1[, l] ~ g2))$adj.r.squared
                if(!is.na(r2_)) {
                    r2_2 = r2_2 + r2_
                    r2_n2 = r2_n2 + 1
                }
            }
            if(r2_n1 > 0 && r2_n2 > 0) {
                r2 = ((r2_1 / r2_n1) + (r2_2 / r2_n2)) / 2
            } else {
                r2 = NaN
            }
        } else {
            r2 = NaN
        }
        result[i-first+1, j] = r2
        k = k + 1
        if(j >= first && j <= last) {
            result[j-first+1, i] = r2
            k = k + 1
        }
        if(k %% 10000 == 0) {
            cat(k, 'genes out of', nrow(result)*ncol(result), '\n', file=stderr())
        }
    }
}

write.table(result, file=stdout(), col.names=(chunk == 1),
            row.names=TRUE, quote=FALSE, sep='\t')
