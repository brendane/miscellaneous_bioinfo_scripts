#!/usr/bin/env Rscript
#
# Conduct an association analysis using QTCAT. Requires a tab-delimited
# file with strain names in the first column and phenotype values in the
# second column (with NA string recognized by R for missing data). Also
# requires a file with genotypes in QTCAT format.
#
# qtcat_gwas.r --alpha <FDR level> [--maf <min MAF> --missing <max missing>]
#   --family <gaussian or binomial> <pheno> <snps>
#
# This script was written using the QTCAT documentation on github as a
# guide (https://github.com/QTCAT/qtcat.data).
#

library(optparse)

library(qtcat)

optlist = list(make_option('--output'),
               make_option('--nsplits', type='double', default=50),
               make_option('--alpha', type='double', default=NaN),
               make_option('--maf', type='double', default=0),
               make_option('--family'),
               make_option('--keep'),
               make_option('--missing', type='double', default=1.0))
parser = OptionParser(option_list=optlist)
opts = parse_args(parser, positional_arguments=TRUE)

outfile = opts$options$output
n_splits = opts$options$nsplits
alpha = opts$options$alpha
min_maf = opts$options$maf
max_miss = opts$options$missing
model_family = opts$options$family
keep_file = opts$options$keep
pheno_file = opts$args[1]
var_file = opts$args[2]


### Read data
var_data = read.snpData(var_file, sep='\t', quote='', na.string='NN')
pheno_data = read.table(pheno_file, sep='\t', as.is=TRUE,
                        colClasses=c('character', 'numeric'))

if(!is.null(keep_file)) {
    keep = scan(keep_file, what='character')
    var_data = var_data[rownames(var_data) %in% keep, ]
    pheno_data = pheno_data[pheno_data[, 1] %in% keep, ]
}

non_missing_strains = pheno_data[!is.na(pheno_data[, 2]), 1]
var_data = var_data[rownames(var_data) %in% non_missing_strains, ]
pheno_data = pheno_data[!is.na(pheno_data[, 2]), ]


### Filter variants
missing_freq = naFreq(var_data, 2)
maf = alleleFreq(var_data, maf=TRUE)
filt_var_data = var_data[, missing_freq <= max_miss & maf >= min_maf]

### Counts of data left after filtering
cat('Number of variants =', ncol(filt_var_data), '\n', file=stderr())
cat('Number of individuals =', nrow(filt_var_data), '\n', file=stderr())
cat('Number of individuals phenotyped =', nrow(pheno_data), '\n', file=stderr())
category_counts = table(pheno_data[, 2])   # Really only applies to binary traits
cat('Number of categories =', length(category_counts), '\n', file=stderr())
cat('Number in categories =', category_counts, '\n', file=stderr())


### Cluster variants
var_clust = qtcatClust(snp=filt_var_data)


### Create genotype and phenotype objects
pheno_obj = qtcatPheno(names=pheno_data[, 1],
                       pheno=pheno_data[, 2],
                       family=model_family)
geno_obj = qtcatGeno(snp=filt_var_data,
                     snpClust=var_clust)


### Run "HIT" analysis
hit_result = qtcatHit(pheno_obj, geno_obj, B=n_splits)
qtcs = qtcatQtc(hit_result, alpha=alpha)


### Write output
write.table(qtcs, file=outfile, sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
