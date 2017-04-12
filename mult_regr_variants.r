#!/usr/bin/env Rscript
#
# This script does multiple linear regression stuff with sets of variants
# and phenotypes.
#

library(magrittr)
library(data.table)
library(dplyr)
library(dtplyr)
library(MuMIn)

impute = function(x) {
    apply(x, 2, function(x) {
          m = mean(x, na.rm=TRUE)
          x[is.na(x)] = m
          x
})
}


#==============================================================================#

cargs = commandArgs(trailingOnly=TRUE)
geno_data_file = cargs[1] # .tsv file created from VCF with 0 and 1 for genotypes
                          # can have or not have an "rs" column
top_vars_file = cargs[2]  # File with rs; rows are replicates, columns are ordered
pheno_file = cargs[3]     # File with two columns: strain and phenotype value, no
                          # header
geno_sampling = as.numeric(cargs[4]) # Number of variants to sample to get genomic pve
calc_genomic_pve = as.logical(as.numeric(cargs[5]))  # Whether or not to calculate the genomic pve
do_interaction = as.logical(as.numeric(cargs[6]))  # Whether or not to run the interaction code
output_prefix = cargs[7]

## Get the genotype data
## It seems easiest to grab this from the tsv file that was used as
## input for the LD script rather than trying to use the VCF or plink
## files.
geno_data = fread(geno_data_file, header=TRUE)
if(!('rs' %in% colnames(geno_data))) {
    geno_data = geno_data %>% mutate(rs=paste0(contig, '-', pos))
}
geno_data %<>% setkey(cols='rs')

## Get the phenotype data and put it in the same order as the
## genotype data
pheno = read.table(pheno_file, as.is=TRUE)
#if(!all(pheno[, 1] %in% colnames(geno_data))) {
#    print(pheno[!pheno[, 1] %in% colnames(geno_data), 1])
#    stop('Some strains are missing from the the genotype data')
#}
pheno = pheno[match(colnames(geno_data), pheno[, 1]), ]
pheno = pheno[!is.na(pheno[, 1]), ]

geno_data = geno_data[, c('contig', 'pos', 'rs', pheno[, 1]), with=FALSE]
geno_mat = geno_data %>% select(-contig, -pos, -rs) %>%
    as.matrix()

## Read in list of top variants
top_vars = read.table(top_vars_file, as.is=TRUE)

if(calc_genomic_pve) {
    sampled_vars = sample(1:nrow(geno_data), geno_sampling, FALSE)
}

## Loop over replicates
for(i in 1:nrow(top_vars)) {

    ## List of variants for this replicate
    focal_vars = as.character(top_vars[i, ])
    n_vars = length(focal_vars)

    ## Get the genotypes for this replicate and impute missing data by
    ## making missing data equal the mean genotype; this is the
    ## strategy used by GEMMA, and it avoids issues with losing lots of
    ## individuals when there is missing data without biasing the results
    ## much.
    top_vars_geno = geno_data[focal_vars, on='rs'] %>%
        select(-contig, -pos, -rs) %>%
        as.matrix %>%
        t() %>%
        impute()

    ## Each variant by itself and each variant added to the model
    ## cumulatively
    output_each = data.frame('repl'=rep(i, n_vars),
                             'var'=focal_vars,
                             'r2'=numeric(n_vars)*NaN,
                             'aic'=numeric(n_vars)*NaN,
                             'r2adj'=numeric(n_vars)*NaN,
                             'geno_pve'=numeric(n_vars)*NaN)
    output_cuml = data.frame('repl'=rep(i, n_vars),
                             'var'=focal_vars,
                             'r2'=numeric(n_vars)*NaN,
                             'aic'=numeric(n_vars)*NaN,
                             'r2adj'=numeric(n_vars)*NaN,
                             'geno_pve'=numeric(n_vars)*NaN)

    ## Loop over variants
    cuml_cols = logical(n_vars)
    for(j in 1:n_vars) {
        ## If there is no non-missing data, continue
        if(all(is.na(top_vars_geno[, j]))) {
            cuml_cols[j] = FALSE
            next
        }
        ## Do the regressions
        cuml_cols[j] = TRUE
        m_cuml = lm(pheno[, 2] ~ top_vars_geno[, cuml_cols])
        sum_cuml = summary(m_cuml)
        m_each = lm(pheno[, 2] ~ top_vars_geno[, j])
        sum_each = summary(m_each)
        output_each$r2[j] = sum_each$r.squared
        output_each$r2adj[j] = sum_each$adj.r.squared
        output_each$aic = AICc(m_each)
        output_cuml$r2[j] = sum_cuml$r.squared
        output_cuml$r2adj[j] = sum_cuml$adj.r.squared
        output_cuml$aic = AICc(m_cuml)

        ## How much genomic variation is predicted by these variants
        if(calc_genomic_pve) {
            r2_genome_sum_each = 0
            r2_genome_count_each = 0
            r2_genome_sum_cuml = 0
            r2_genome_count_cuml = 0
            for(k in seq_along(sampled_vars)) {
                r2_each_ = summary(lm(as.numeric(geno_mat[sampled_vars[k], ]) ~ top_vars_geno[, j]))$r.squared
                r2_cuml_ = summary(lm(as.numeric(geno_mat[sampled_vars[k], ]) ~ top_vars_geno[, cuml_cols]))$r.squared
                if(!is.na(r2_each_)) {
                    r2_genome_sum_each %<>% + r2_each_
                    r2_genome_count_each %<>% + 1
                }
                if(!is.na(r2_cuml_)) {
                    r2_genome_sum_cuml %<>% + r2_cuml_
                    r2_genome_count_cuml %<>% + 1
                }
            }
            output_each$geno_pve[j] = r2_genome_sum_each / r2_genome_count_each
            output_cuml$geno_pve[j] = r2_genome_sum_cuml / r2_genome_count_cuml
        }
    }

    write.table(output_cuml,
                file=paste0(output_prefix, '.cumulative.tsv'),
                col.names=(i == 1), row.names=FALSE, quote=FALSE, sep='\t',
                append=(i != 1))
    write.table(output_each,
                file=paste0(output_prefix, '.each.tsv'),
                col.names=(i == 1), row.names=FALSE, quote=FALSE, sep='\t',
                append=(i != 1))

    if(do_interaction) {

        ## Pairwise output
        n_pairs = choose(n_vars, 2)
        output_pair = data.frame('repl'=rep(i, n_pairs),
                                 'var1'=character(n_pairs),
                                 'var2'=character(n_pairs),
                                 'r2_additive'=numeric(n_pairs),
                                 'r2_interact'=numeric(n_pairs),
                                 'r2adj_additive'=numeric(n_pairs),
                                 'r2adj_interact'=numeric(n_pairs),
                                 'aic_additive'=numeric(n_pairs),
                                 'aic_interact'=numeric(n_pairs),
                                 'LRT_p'=numeric(n_pairs),
                                 stringsAsFactors=FALSE)

        ## Loop over pairs, testing for interaction
        ## Use AICc function because the sample size is a little small
        pair_list = combn(1:n_vars, 2)
        for(j in 1:ncol(pair_list)) {
            a = pair_list[1, j]
            b = pair_list[2, j]
            var1 = focal_vars[a]
            var2 = focal_vars[b]
            m0 = lm(pheno[, 2] ~ top_vars_geno[, a] + top_vars_geno[, b])
            m1 = lm(pheno[, 2] ~ top_vars_geno[, a] * top_vars_geno[, b])
            s0 = summary(m0)
            s1 = summary(m1)
            output_pair[j, 'var1'] = var1
            output_pair[j, 'var2'] = var2
            output_pair[j, 'r2_additive'] = s0$r.squared
            output_pair[j, 'r2adj_additive'] = s0$adj.r.squared
            output_pair[j, 'r2_interact'] = s1$r.squared
            output_pair[j, 'r2adj_interact'] = s1$adj.r.squared
            output_pair[j, 'aic_additive'] = AICc(m0)
            output_pair[j, 'aic_interact'] = AICc(m1)
            output_pair[j, 'LRT_p'] = anova(m0, m1, test='LRT')[['Pr(>Chi)']][2]
        }
        write.table(output_pair,
                    file=paste0(output_prefix, '.interaction.tsv'),
                    col.names=(i == 1), row.names=FALSE, quote=FALSE, sep='\t',
                    append=(i != 1))
    }

}
