#!/usr/bin/env Rscript
#
# Conduct an association analysis using QTCAT. Requires a tab-delimited
# file with strain names in the first column and phenotype values in the
# second column (with NA string recognized by R for missing data). Also
# requires a file with genotypes in QTCAT format.
#
# qtcat_gwas.r --alpha <FDR level> [--maf <min MAF> --missing <max missing>]
#   <pheno> <snps>
#
# This script was written using the QTCAT documentation on github as a
# guide (https://github.com/QTCAT/qtcat.data).
#

library(optparse)

library(qtcat)

optlist = list(make_option('--output'),
               make_option('--alpha', type='double', default=NaN),
               make_option('--maf', type='double', default=0),
               make_option('--missing', type='double', default=1.0))
parser = OptionParser(option_list=optlist)
opts = parse_args(parser, positional_arguments=TRUE)

alpha = opts$options$alpha
min_maf = opts$options$maf
max_miss = opts$options$missing
pheno_file = opts$args[1]
var_file = opts$args[2]


### Read data
var_data = read.snpData(var_file)


### Filter variants


### Cluster variants


### Create genotype and phenotype objects


### Run "HIT" analysis
