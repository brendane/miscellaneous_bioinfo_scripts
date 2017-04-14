#!/usr/bin/env Rscript
#
# Use generalized linear models to test for differences in allele
# frequency.
#
# Each read is treated as an observation, and pools are the replicates.
# A likelihood ratio test is used to produce a p-value.
#
# Usage:
#
#   snp_by_snp_glm.r [--randomize] <sync file> <file with pool names>
#       <group 1 pools> <group 2 pools>
#
#   The names of the pools in each group should be separated by commas.
#
#==============================================================================#


library(data.table)

options(warn=1)

#==============================================================================#

## From a sync file row, get read counts.
get_counts = function(sync_row, pools, pool_names) {
    contig = sync_row[[1]]
    pos = sync_row[[2]]
    ref_allele = toupper(sync_row[[3]])
    counts = structure(vector('list', length=length(pools)),
                       names=pools)
    for(pn in pools) {
        x = sync_row[[which(pool_names == pn) + 3]]
        cnts = as.numeric(unlist(strsplit(x, ':')))
        names(cnts) = c('A', 'T', 'C', 'G', 'N', 'D')
        ref_cnt = cnts[ref_allele]
        alt_cnt = (sum(cnts) - cnts['N']) - ref_cnt
        counts[[pn]] = c('ref'=ref_cnt, 'alt'=alt_cnt)
    }
    counts
}

## Test if this set of counts is variable
is_variable = function(counts) {
    ref = sum(sapply(counts, function(x) x[1]))
    alt = sum(sapply(counts, function(x) x[2]))
    return((ref != 0) && (alt != 0))
}

null_result = list('p'=NaN, 'effect'=NaN, term=NA, converged=NA,
                   'freqs0'=NA, 'freqs1'=NA)

## Run a GLM Use a likelihood ratio test to get p-values.
## Counts needs to be in the same order as groups
test = function(counts, groups) {
    resp = matrix(unlist(counts), ncol=2, byrow=TRUE)
    pools = factor(names(counts))
    trt = rep(names(groups), sapply(groups, length))
    freqs = structure(round(sapply(counts, function(x) x[1] / sum(x)), 4),
                      names=names(counts))
    freqs0 = paste(freqs[names(freqs) %in% groups[[1]]], collapse=',')
    freqs1 = paste(freqs[names(freqs) %in% groups[[2]]], collapse=',')
    ret = null_result
    if(length(unique(trt)) > 1) {
        m = glm(resp ~ trt, family=binomial('logit'))
        mnull = glm(resp ~ 1, family=binomial('logit'))
        lrt_p = anova(m, mnull, test='Chisq')[['Pr(>Chi)']][2]
        s = summary(m)
        conv = (m$converged && mnull$converged)
        fe = coefficients(m)[2]
        ret = list('p'=lrt_p,'effect'=fe, 'term'=names(fe),
                   'converged'=as.numeric(conv),
                   'freqs0'=freqs0, 'freqs1'=freqs1)
    }
    ret
}

#==============================================================================#

## Get arguments to the script
cargs = commandArgs(trailingOnly=TRUE)
i = 0
rand = FALSE
if(cargs[1] == '--randomize') {
    rand = TRUE
    i = 1
}
sync_file_name = cargs[i + 1]
pool_file_name = cargs[i + 2]
group_1 = Filter(function(x) x != '', unlist(strsplit(cargs[i + 3], ',')))
group_2 = Filter(function(x) x != '', unlist(strsplit(cargs[i + 4], ',')))

## Read in data
read_counts = fread(sync_file_name, sep='\t', stringsAsFactors=FALSE,
                    verbose=FALSE, showProgress=FALSE)
pool_names = scan(pool_file_name, what='character')

if(rand) {
    all_pools = c(group_1, group_2)
    rand_group_1 = sample(all_pools, length(group_1), FALSE)
    rand_group_2 = all_pools[!all_pools %in% rand_group_1]
    trts = list('group1'=rand_group_1, 'group2'=rand_group_2)
} else {
    trts = list('group1'=group_1, 'group2'=group_2)
}

## Test every SNP and print results
cat('contig\tpos\tref\tp\teffect\tterm\tconverged\tfreqs0\tfreqs1\n', file=stdout())
for(i in 1:nrow(read_counts)) {
    row_data = read_counts[i]
    contig = row_data[[1]]
    pos = row_data[[2]]
    ref = row_data[[3]]
    counts = get_counts(row_data, unlist(trts), pool_names)
    if(is_variable(counts)) {
        test_results = tryCatch(test(counts, trts),
                                error=function(e) null_result)
    } else {
        test_results = null_result
    }

    cat(paste(contig, pos, ref, test_results$p,
              test_results$effect, test_results$term, test_results$converged,
              test_results$freqs0, test_results$freqs1,
              sep='\t'),
        '\n', sep='', file=stdout())
}
