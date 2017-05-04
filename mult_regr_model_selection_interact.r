#!/usr/bin/env Rscript
#
# This script does simple forward model selection on variants and
# phenotypes. A simple imputation method (just taking the mean) is used
# to deal with missing genotype data.
#

library(magrittr)
library(data.table)
library(dplyr)
library(dtplyr)
library(MuMIn)

## Function to impute missing data by replacing missing genotypes with
## the mean. Expects alleles to be coded as 0 and 1.
impute = function(x) {
    apply(x, 2, function(x) {
          m = mean(x, na.rm=TRUE)
          x[is.na(x)] = m
          x
})
}

#==============================================================================#

cargs = commandArgs(trailingOnly=TRUE)
if(length(cargs) < 3) {
    cat('mult_regr_model_selections.r <genotypes> <variants> <phenotypes>\n\n',
        file=stderr())
}
geno_data_file = cargs[1] # .tsv file created from VCF with 0 and 1 for genotypes
                          # can have or not have an "rs" column
top_vars_file = cargs[2]  # File with rs; rows are replicates, columns are ordered
pheno_file = cargs[3]     # File with two columns: strain and phenotype value, no
                          # header

## Get the genotype data
## It seems easiest to grab this from the tsv file that was used as
## input for the LD script rather than trying to use the VCF or plink
## files. If there isn't an "rs" column (variant identification), it
## joins the chromosome and position with "-" to get an ID.
geno_data = fread(geno_data_file, header=TRUE,
                  verbose=FALSE, showProgress=FALSE)
if(!('rs' %in% colnames(geno_data))) {
    geno_data = geno_data %>% mutate(rs=paste0(contig, '-', pos))
}
geno_data %<>% setkey(cols='rs')

## Get the phenotype data and put it in the same order as the
## genotype data
pheno = read.table(pheno_file, as.is=TRUE)
pheno = pheno[match(colnames(geno_data), pheno[, 1]), ]
pheno = pheno[!is.na(pheno[, 1]), ]

## Make a matrix with just the genotypes
geno_data = geno_data[, c('contig', 'pos', 'rs', pheno[, 1]), with=FALSE]
geno_mat = geno_data %>% select(-contig, -pos, -rs) %>%
    as.matrix()

## Read in list of top variants
top_vars = read.table(top_vars_file, as.is=TRUE)

## Loop over replicates
cat('replicate\tnvars\tR2\tAIC\tbest\tadded\tvars\n')
for(i in 1:nrow(top_vars)) {

    ## List of variants for this replicate
    focal_vars = as.character(top_vars[i, ])

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
    all_missing = logical(length(focal_vars))
    for(j in 1:length(focal_vars)) {
        if(all(is.na(top_vars_geno[, j]))) {
            all_missing[j] = TRUE
        }
    }
    top_vars_geno = top_vars_geno[, !all_missing]
    focal_vars = focal_vars[!all_missing]
    focal_interacts = t(combn(focal_vars, 2))
    n_vars = length(focal_vars)
    n_interacts = length(focal_interacts)
    n_terms = n_vars + n_interacts

    ## Each variant by itself
    indiv_vars = data.frame('var'=focal_vars,
                            'r2'=numeric(n_vars),
                            'aic'=numeric(n_vars),
                            stringsAsFactors=FALSE)
    for(j in 1:n_vars) {
        m = lm(pheno[, 2] ~ top_vars_geno[, j])
        indiv_vars[j, 'r2'] = summary(m)$r.squared
        indiv_vars[j, 'aic'] = AICc(m)
        indiv_vars[j, 'var'] = focal_vars[j]
    }
    write.table(data.frame('replicate'=rep(i, n_vars),
                           'nvars'=rep(1, n_vars),
                           'R2'=indiv_vars['r2'],
                           'AIC'=indiv_vars['aic'],
                           'best'=rep(NA, n_vars),
                           'added'=indiv_vars['var'],
                           'vars'=indiv_vars['var']),
                file=stdout(), col.names=FALSE, sep='\t',
                row.names=FALSE, quote=FALSE)

    ## Model selection
    sel_methods = c('R2', 'AIC')
    n_meth = length(sel_methods)
    model_selection = data.frame(replicate=rep(i, n_terms*n_meth),
                                 'nvars'=rep(1:n_terms, each=n_meth),
                                 'r2'=numeric(n_terms*n_meth),
                                 'aic'=numeric(n_terms*n_meth),
                                 'best'=rep(1, n_terms*n_meth),
                                 'added'=character(n_terms*n_meth),
                                 'vars'=character(n_terms*n_meth),
                                 stringsAsFactors=FALSE)
    for(sm in seq_along(sel_methods)) {
        method = sel_methods[sm]
        model_vars = logical(n_vars)
        model_i_vars = logical(n_interacts)
        for(j in 1:n_vars) {
            tested_vars = which(!model_vars)
            if(sum(model_vars) > 1) {
                ## TODO: interaction terms to test need to be variables
                ## currently in the model as main effects
                tested_interactions = which(!model_i_vars & ...)
            }
            ## TODO - needs to include interactions
            models = vector('list', length=sum(!model_vars))
            for(k in seq_along(tested_vars)) {
                cols = model_vars
                cols[tested_vars[k]] = TRUE
                models[[k]] = lm(pheno[, 2] ~ top_vars_geno[, cols])
            }
            ## This part is mostly okay
            r2s = sapply(models, function(m) summary(m)$r.squared)
            aics = sapply(models, function(m) AICc(m))
            if(method == 'R2') {
                best = which(r2s == max(r2s))
                best = which(r2s == max(r2s))[1]
            } else if(method == 'AIC') {
                best = which(aics == min(aics))
                best = which(aics == min(aics))[1]
            }
            if(length(best) > 1) {
                ## This shouldn't happen if the data has been thinned
                ## by LD group properly
                cat('Warning: more than one best choice\n', file=stderr())
                best = best[1]
            }
            model_vars[tested_vars[best]] = TRUE
            ## TODO adjust index n_vars*... for interactions
            model_selection[(n_vars*(sm-1))+j, 'nvars'] = j
            model_selection[(n_vars*(sm-1))+j, 'r2'] = r2s[best]
            model_selection[(n_vars*(sm-1))+j, 'aic'] = aics[best]
            model_selection[(n_vars*(sm-1))+j, 'aic'] = method
            model_selection[(n_vars*(sm-1))+j, 'added'] = focal_vars[tested_vars[best]]
            ## TODO - this line needs to handle interactions
            model_selection[(n_vars*(sm-1))+j, 'vars'] = paste(focal_vars[model_vars],
                                                               collapse=',')
        }
    }
    write.table(model_selection, file=stdout(), sep='\t',
                col.names=FALSE, row.names=FALSE, quote=FALSE)
}
