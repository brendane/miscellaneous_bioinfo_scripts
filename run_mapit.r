#!/usr/bin/env Rscript
#
# Script to run MAPIT -- an R program that calculates marginal epistasis
# values for variants.
#
# run_mapit.r <directory with MAPIT code> <genotypes> <phenotypes>
#
# genotypes is assumed to be tab-delimited with chrom, pos, and rs columns
# before the genotypes. The genotypes should be 1 or 0, and the file should
# have column names.
#
# phenotypes is a two column file, with no header. The first column should
# be the strain and the second column should be the phenotype value.
#

library(data.table)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)

cargs = commandArgs(trailingOnly=TRUE)
mapit_dir = cargs[1]
geno_file = cargs[2]
pheno_file = cargs[3]

min_maf = 0.05

source(file.path(mapit_dir, 'MAPIT.R'))
sourceCpp(file.path(mapit_dir, 'MAPIT.cpp'))
source(file.path(mapit_dir, 'QQPlot.R'))

## Genotypes and minor allele frequency filtering
genos = fread(geno_file, header=TRUE)
X = as.matrix(genos[, -(1:3)])
rownames(X) = genos[['rs']]

maf = apply(X, 1, function(x) {
            af = sum(x, na.rm=TRUE) / sum(!is.na(x))
            if(af > 0.5) {
                1-af
            } else {
                af
            }
})

X = X[maf >= min_maf, ]

## Phenotypes
phenos = read.csv(pheno_file, sep='\t', header=FALSE, as.is=TRUE,
                  colClasses=c('character', 'numeric'))
Y = phenos[as.character(colnames(X)), 2]


mapit = MAPIT(X, Y)
hybrid_pvals = mapit[['pvalues']]
names(hybrid_pvals) = rownames(X)

