#!/usr/bin/env Rscript
#
# Script to run MAPIT -- an R program that calculates marginal epistasis
# values for variants.
#
# run_mapit.r <directory with MAPIT code> <genotypes> <phenotypes> <output prefix>
#
# genotypes is assumed to be tab-delimited with chrom, pos, and rs columns
# before the genotypes. The genotypes should be 1 or 0, and the file should
# have column names. Missing data is imputed by taking the average.
#
# phenotypes is a two column file, with no header. The first column should
# be the strain and the second column should be the phenotype value.
#
# The p-values are calculated with the Davies method, which is the most
# accurate, though slowest, method available in MAPIT.
#

library(data.table)
library(CompQuadForm)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(optparse)

optlist = list(make_option('--output'),
               make_option('--maf', type='double', default=0.05),
               make_option('--extract'))
opts = parse_args(OptionParser(option_list=optlist),
                  positional_arguments=TRUE)
mapit_dir = opts[['args']][1]
geno_file = opts[['args']][2]
pheno_file = opts[['args']][3]
outprefix = opts[['options']][['output']]
min_maf = opts[['options']][['maf']]
extract_file = opts[['options']][['extract']]

source(file.path(mapit_dir, 'MAPIT.R'))
sourceCpp(file.path(mapit_dir, 'MAPIT.cpp'))
source(file.path(mapit_dir, 'QQPlot.R'))


## Genotypes and minor allele frequency filtering
genos = fread(geno_file, header=TRUE)
X = as.matrix(genos[, -(1:3)])
rownames(X) = genos[['rs']]

if(!is.null(extract_file)) {
    var_list = scan(extract_file, what='character')
    var_list = var_list[var_list %in% rownames(X)]
    X = X[var_list, ]
}

maf = apply(X, 1, function(x) {
            af = sum(x, na.rm=TRUE) / sum(!is.na(x))
            if(af > 0.5) {
                1-af
            } else {
                af
            }
})

X = X[maf >= min_maf, ]

X_imp = t(apply(X, 1, function(x) {
                m = mean(x, na.rm=TRUE)
                ifelse(is.na(x), m, x)
}))

## Phenotypes
phenos = read.csv(pheno_file, sep='\t', header=FALSE, as.is=TRUE,
                  colClasses=c('character', 'numeric'))
rownames(phenos) = phenos[, 1]
Y = phenos[as.character(colnames(X_imp)), 2]

## Get rid of individuals with missing data
X_filt = X_imp[, !is.na(Y)]
Y_filt = Y[!is.na(Y)]

cat('Running MAPIT\n', file=stderr())

## Run mapit with Davie's method for calculating p-values
mapit = MAPIT(X_filt, Y_filt, hybrid=FALSE, test='davies')

## Write output files
png(paste0(outprefix,'.qqplots.png'))
ggd.qqplot(mapit[['pvalues']])

write.table(data.frame('rs'=rownames(X_filt),
                       'p'=mapit[['pvalues']],
                       'pve'=mapit[['pves']]),
            file=paste0(outprefix, '.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
