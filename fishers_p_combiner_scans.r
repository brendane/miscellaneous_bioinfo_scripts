#!/usr/bin/env Rscript
#
# Combine results from genome scans using Fisher's method.
#
# Assumes the input files are tab-delimited with no headers and that the
# 1st column is the chromosome, the 2nd column is the position, the
# 3rd column is the p-value, and the 4th column is whether or not the
# site is a fixed difference. Additional columns are fine.
#

library(data.table)

fisher_combine_p = function(ps) {
        chi_sq = -2 * sum(log(ps[!is.na(ps)]))
    df = 2 * sum(!is.na(ps))
        pchisq(chi_sq, df, lower.tail=FALSE)
}

check_site_matches = function(cmat, pmat) {
    pass = TRUE
    pass = pass && all(apply(cmat, 1,
                             function(x) all(x[1] == x[2:length(x)])))
    pass = pass && all(apply(pmat, 1,
                             function(x) all(x[1] == x[2:length(x)])))
    pass
}

cargs = commandArgs(trailingOnly=TRUE)

data_list = lapply(lapply(cargs, fread, sep='\t'),
                   function(x) { names(x)[1:4] = c('chrom', 'pos', 'p', 'fix'); x })
f_matrix = sapply(data_list, function(x) x[['fix']])
p_matrix = sapply(data_list, function(x) x[['p']])
colnames(p_matrix) = basename(dirname(cargs))
chr_matrix = sapply(data_list, function(x) x[['chrom']])
pos_matrix = sapply(data_list, function(x) x[['pos']])

if(!check_site_matches(chr_matrix, pos_matrix))
    stop('Sites and positions do not match')

combined_ps = apply(p_matrix, 1, fisher_combine_p)
combined_fs = rowSums(f_matrix, na.rm=TRUE)

write.table(data.frame(chrom=chr_matrix[, 1], pos=pos_matrix[, 1],
                       p=combined_ps, fixed=combined_fs,
                       p_matrix),
            file=stdout(), sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
