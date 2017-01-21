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

library(data.table)

#==============================================================================#

cargs = commandArgs(trailingOnly=TRUE)

pool_freqs_file = cargs[1]
base_calls_file = cargs[2]

pool_freqs = fread(pool_freqs_file, sep='\t', header=TRUE)
base_calls = fread(base_calls_file, sep='\t', header=TRUE)
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

# alternate allele frequency as response, strain base calls as
# predictors, no intercept.
model_formula = paste('ref_freqs ~ (-1) + ',
                      paste(strains, collapse=' + '))

# Linear regression:
model_results = lm(as.formula(model_formula), data=base_calls,
                   weights=pool_freqs[['n_reads']], na.action=na.exclude)
sm = summary(model_results)
output_table = sapply(strains, function(x) sm$coefficients[x, 1:2])
output_table = data.frame(stat = rownames(output_table),
                          output_table,
                          R2 = c(sm$adj.r.squared, ''))
write.table(output_table, file=stdout(), col.names=TRUE,
            row.names=FALSE, quote=FALSE, sep='\t')

# GLM
model_results = glm(as.formula(model_formula), data=base_calls,
                    family=binomial('logit'), weights=pool_freqs[['n_reads']],
                    na.action=na.exclude)
# And try to adjust the coefficients back to original scale
cf = exp(model_results$coefficients) / sum(exp(model_results$coefficients))
output_table = data.frame(stat = 'Estimate-GLM', cf, R2 = '')
write.table(output_table, file=stdout(), col.names=FALSE,
            row.names=FALSE, quote=FALSE, sep='\t')

# Optimization method; minimize sums of squares with no transformation,
# and constrain estimates to between 0 and 1. Another faster way
# to accomplish the constraint is to modify get_ss to strongly penalize
# estimates that are not in the right range, but did not perform as
# well in testing with fake data.
get_ss = function(par, calls, freqs) {
    ## cfs = coefficients (estimated strain frequencies)
    ## calls = matrix: columns are strains, rows are sites; 1 or 0
    ## freqs = vector of pool allele frequencies (one element per site)
    ## Value: sum of squares
    prd = (calls %*% par)[, 1]
    ss = sum((prd - freqs)^2)
    ss
}

bc = as.matrix(base_calls)[, -(1:2)]
storage.mode(bc) = 'numeric'
nmiss = rowSums(is.na(bc)) 
bc = bc[nmiss == 0 & !is.na(ref_freqs), ]
rf = ref_freqs[nmiss == 0 & !is.na(ref_freqs)]
o = optim(par=rep(1/ncol(bc), ncol(bc)), get_ss, calls=bc,
          freqs=rf, lower=0, upper=1, method='L-BFGS-B')
write.table(data.frame(stat='Estimate', o$par, 'R2'=''),
            file=stdout(), row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep='\t')
