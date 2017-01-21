#!/usr/bin/env Rscript
#
# Use generalized linear mixed models to test for differences in allele
# frequency.
#
# Each read is treated as an observation, and pools are treated as random
# effects. A likelihood ratio test is used to produce a p-value.
#
# Usage:
#
#   snp_by_snp_glmm_test.r <sync file> <file with pool names>
#       <group 1 pools> <group 2 pools>
#
#   The names of the pools in each group should be separated by commas.
#
#==============================================================================#


library(data.table)
library(lme4)

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
                   'wald_p'=NaN)

## Run a GLMM with pool as a random effect. Use the default
## Wald Z-test and likelihood ratio test to get p-values
test = function(counts, groups) {
    resp = unlist(sapply(counts, function(x) c(rep(0, x[1]), rep(1, x[2]))))
    n_reads = sapply(counts, sum)
    pools = factor(rep(names(counts), n_reads))
    trt = rep(names(groups),
              sapply(groups, function(x) sum(n_reads[x])))
    ret = null_result
    if(length(unique(trt)) > 1) {
        ## Random intercept for each pool
        m = glmer(resp ~ trt + (1|pools), family=binomial('logit'))
        mnull = glmer(resp ~ (1|pools), family=binomial('logit'))
        lrt_p = anova(m, mnull)[['Pr(>Chisq)']][2]
        s = summary(m)
        conv = (length(s$optinfo$conv$lme4) == 0)
        p = s$coefficients[2, 4]
        fe = fixef(m)[2]
        ret = list('p'=lrt_p, 'wald_p'=p, 'effect'=fe, 'term'=names(fe),
                   'converged'=as.numeric(conv))
    }
    ret
}

#==============================================================================#

## Get arguments to the script
cargs = commandArgs(trailingOnly=TRUE)
sync_file_name = cargs[1]
pool_file_name = cargs[2]
group_1 = Filter(function(x) x != '', unlist(strsplit(cargs[3], ',')))
group_2 = Filter(function(x) x != '', unlist(strsplit(cargs[4], ',')))

## Read in data
read_counts = fread(sync_file_name, sep='\t', stringsAsFactors=FALSE)
pool_names = scan(pool_file_name, what='character')

trts = list('group1'=group_1, 'group2'=group_2)

## Test every SNP and print results
cat('contig\tpos\tref\tp\twald_p\teffect\tterm\tconverged\n', file=stdout())
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
    cat(paste(contig, pos, ref, test_results$p, test_results$wald_p,
              test_results$effect, test_results$term, test_results$converged,
              sep='\t'),
        '\n', sep='', file=stdout())
}
