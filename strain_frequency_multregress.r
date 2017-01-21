#!/usr/bin/env Rscript
#
# Use multiple regression to get strain frequencies.
#
# strain_frequency_multregress.r <pool freqs> <base calls>
#
# Note that the base calls file should have a 1 for reference and 0 for
# alternate allele call, which is the opposite of a VCF file.
#
#==============================================================================#

suppressMessages(library(data.table))

#==============================================================================#

get_ss = function(par, calls, freqs) {
    ## cfs = coefficients (estimated strain frequencies)
    ## calls = matrix: columns are strains, rows are sites; 1 or 0
    ## freqs = vector of pool allele frequencies (one element per site)
    ## Value: sum of squares
    prd = (calls %*% par)[, 1]
    ss = sum((prd - freqs)^2)
    ss
}

get_sse = function(par, calls, freqs) {
        ## cfs = coefficients (estimated strain frequencies)
        ## calls = matrix: columns are strains, rows are sites; 1 or 0
        ## freqs = vector of pool allele frequencies (one element per site)
        ## Value: sum of squares
        prd = (calls %*% par)[, 1]
    ss = sum((prd - freqs)^2)
        ss
}

get_fitted = function(par, calls) {
        (calls %*% par)[, 1]
}

get_sst = function(freqs) {
        sum((mean(freqs) - freqs)^2)
}

get_r2 = function(par, calls, freqs) {
        ## par = $par slot from optimization - the coefficients
        ## calls = a matrix of allele calls (columns = strains)
        ## freqs = the pool allele frequencies
        sse = get_sse(par, calls, freqs)
    sst = get_sst(freqs)
        n = nrow(calls)
        m = length(par)
            v = n - m
            1 - ( (sse * (n - 1)) / (sst * v) )
}

#==============================================================================#

cargs = commandArgs(trailingOnly=TRUE)

pool_freqs_file = cargs[1]
base_calls_file = cargs[2]

sink(stderr())
pool_freqs = suppressMessages(fread(pool_freqs_file, sep='\t',
                                    header=TRUE, verbose=FALSE))
base_calls = suppressMessages(fread(base_calls_file, sep='\t',
                                    header=TRUE, verbose=FALSE))
sink()
names(base_calls)[-(1:2)] = paste0('S', names(base_calls)[-(1:2)])
strains = names(base_calls)[-(1:2)]

base_calls = base_calls[order(base_calls[['contig']], base_calls[['pos']]), ]
pool_freqs = pool_freqs[order(pool_freqs[['contig']], pool_freqs[['pos']]), ]
pool_key = paste0(pool_freqs[['contig']], '-', pool_freqs[['pos']])
bc_key = paste0(base_calls[['contig']], '-', base_calls[['pos']])
base_calls = base_calls[bc_key %in% pool_key, ]
pool_freqs = pool_freqs[pool_key %in% bc_key, ]

if(!all(pool_freqs[['contig']] == base_calls[['contig']]))
    stop('Files are not in the same order')
if(!all(pool_freqs[['pos']] == base_calls[['pos']]))
    stop('Files are not in the same order')


ref_freqs = pool_freqs[['ref_freq']]
ref_freqs[pool_freqs[['n_alleles']] == 0] = NaN

# Reference allele frequency as response, strain base calls as
# predictors, no intercept.
model_formula = paste('ref_freqs ~ (-1) + ',
                      paste(strains, collapse=' + '))

# Linear regression:
cat('Running linear regression\n', file=stderr())
model_results = lm(as.formula(model_formula), data=base_calls,
                   weights=pool_freqs[['n_reads']], na.action=na.exclude)
sm = summary(model_results)
output_table = sapply(strains, function(x) sm$coefficients[x, 1:2])
output_table = data.frame(stat = paste('LM', rownames(output_table), sep='-'),
                          output_table,
                          R2 = c(sm$adj.r.squared, ''))
write.table(output_table, file=stdout(), col.names=TRUE,
            row.names=FALSE, quote=FALSE, sep='\t')

# Optimization method with coefficients constrained between 0 and
# 1
cat('Running constrained optimization\n', file=stderr())
bc = as.matrix(base_calls)[, -(1:2)]
storage.mode(bc) = 'numeric'
nmiss = rowSums(is.na(bc)) 
bc = bc[nmiss == 0 & !is.na(ref_freqs), ]
rf = ref_freqs[nmiss == 0 & !is.na(ref_freqs)]
o = optim(par=rep(1/ncol(bc), ncol(bc)), get_sse, calls=bc,
          freqs=rf, lower=0, upper=1, method='L-BFGS-B')
r2 = get_r2(o$par, bc, rf)
cat('Constrained-Estimate\t', paste(o$par, collapse='\t'), '\t',
    r2, '\n', sep='', file=stdout())
