#!/usr/bin/env Rscript
#
# Script to run MAPIT -- an R program that calculates marginal epistasis
# values for variants.
#
# run_mapit.r <directory with MAPIT code> <genotypes> <phenotypes>
#

library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)

cargs = commandArgs(trailingOnly=TRUE)
mapit_dir = cargs[1]
geno_file = cargs[2]
pheno_file = cargs[3]

source(file.path(mapit_dir, 'MAPIT.R'))
sourceCpp(file.path(mapit_dir, 'MAPIT.cpp'))
source(file.path(mapit_dir, 'QQplot.R'))
