#!/usr/bin/env Rscript
#
# Run Fisher's exact test to compare the allele frequencies in two
# pools.
#
# fishers_exact_pools.r <input1> <input2>
#
# Writes to stdout.
#

#==============================================================================#

cargs = commandArgs(trailingOnly=TRUE)
infile_1 = cargs[1]
infile_2 = cargs[2]

dat_1 = read.csv(infile_1, sep='\t', header=TRUE, as.is=TRUE)
dat_2 = read.csv(infile_2, sep='\t', header=TRUE, as.is=TRUE)

for(i in 1:nrow(dat_1)) {
    d1 = dat_1[i, ]
    d2 = dat_2[i, ]
    if(d1[,'chrom'] != d2[,'chrom'])
        stop(paste('files do not match on line', i))
    if(d1[,'pos'] != d2[,'pos'])
        stop(paste('files do not match on line', i))
    if(d1[,'alleles'] != d2[,'alleles'])
        stop(paste('files do not match on line', i))

    n1 = d1[, 'N']
    n2 = d2[, 'N']

    if(n1 == 0 || n2 == 0) {

        # If there isn't any coverage
        p = NaN
        fixed = NA

    } else {

        # Counts of alleles
        f1 = as.numeric(unlist(strsplit(d1[, 'freqs'], ','))) * n1
        f2 = as.numeric(unlist(strsplit(d2[, 'freqs'], ','))) * n2

        if(length(f1) == 1 || length(f2) == 1) {

            # If there isn't any variation (happens with P/A variants)
            p = NaN
            fixed = NA

        } else {

            # Look for fixed differences
            fixed = 0
            if(any(f1 == n1) && any(f2 == n2) && f2[f1 == n1] == 0)
                fixed = 1

            # Fisher's exact test
            res = fisher.test(rbind(f1, f2), alternative='two.sided')
            p = res$p.value
        }
    }

    cat(d1[, 'chrom'], '\t', d1[, 'pos'], '\t', p, '\t',
        fixed, '\n', sep='', file=stdout())
}
