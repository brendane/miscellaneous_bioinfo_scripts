#!/usr/bin/env Rscript
###
### R script to do run GWAS on pairwise interactions. Choice
### of logistic or gaussian regression, and optional covariates.
###
### DETAILS
###
### For the initial pass, it uses the default settings of lrgpr::glmApply,
### which set missing genotypes to the mean and uses a Wald Chi-square
### test. It also tests all the SNPs - the filters are applied before
### testing interactions.
###
### The second pass uses the regular glm function from R, also with
### default settings, except that NAs are just ignored.
###
#==============================================================================#

library(lrgpr)
library(optparse)

#==============================================================================#

opts = list(make_option('--output'),
            make_option('--test-percent', type='double', default=1),
            make_option('--min-maf', type='double', default=0))
parser = OptionParser(option_list=opts)
argv = parse_args(parser, positional_arguments=TRUE)

outfile = argv$options[['output']]
test_percent = argv$options[['test-percent']]
min_maf = argv$options[['min-maf']]
tfam_file = argv$args[1]
tped_file = argv$args[2]
model_type = argv$args[3]
if(length(argv$args) > 3) {
    covariates_file = argv$args[4]
    covariates = read.table(covariates_file)
    use_covariates = TRUE
} else {
    use_covariates = FALSE
}

## Figure out the type of regression to run
if(model_type == 'gaussian') {
    fam = gaussian(link='identity')
} else if(model_type == 'logistic') {
    fam = binomial(link='logit')
} else {
    stop('Model type not recognized')
}

## Make the binary file
convertToBinary(tped_file, './binary_genotypes', 'TPED')

## Read in data
X = attach.big.matrix('./binary_genotypes_descr', readonly=TRUE)
y = read.table(tfam_file)[, 6]

## Get minor allele frequency and other information
maf = getAlleleFreq(X)
maf[maf > 0.5] = 1 - maf[maf > 0.5]
pos = as.numeric(gsub('.+-', '', colnames(X)))
chr = gsub('-.+', '', colnames(X))

## Set up formula for initial pass
if(use_covariates) {
    f = paste('y ~ ', paste(paste0('covariates[, ', 1:ncol(covariates), ']'),
                            collapse=' + '), '+ SNP')
    f = as.formula(f)
    term = ncol(covariates) + 2
} else {
    f = as.formula('y ~ SNP')
    term = 2
}

## Run initial pass - one SNP at a time
initial_pvalues = glmApply(f, features=X,  terms=term, family=fam)$pValues[, 1]

## Apply filters
snps_to_test = initial_pvalues >= quantile(initial_pvalues,
                                           1-test_percent, na.rm=TRUE)
snps_to_test[is.na(snps_to_test)] = FALSE
snps_to_test = snps_to_test &
                maf >= min_maf
snps_to_test[is.na(snps_to_test)] = FALSE
cat('SNPs to test <=', sum(snps_to_test), '\n',
    file=stderr())

## Run interaction tests
outhandle = file(outfile, 'wb')
cat('chr1\tpos1\tsnp1\tfreq1\t',
    'chr2\tpos2\tsnp2\tfreq2\tconverged\t',
    'N\tp\tbeta\n', sep='', file=outhandle)
for(i in 1:(ncol(X) - 1)) {

    if(!snps_to_test[i])
        next

    chr1 = chr[i]
    pos1 = pos[i]
    maf1 = maf[i]
    snp1 = colnames(X)[i]

    for(j in (i+1):ncol(X)) {

        if(!snps_to_test[j])
            next

        chr2 = chr[j]
        pos2 = pos[j]
        maf2 = maf[j]
        snp2 = colnames(X)[j]

        # Approximate test to avoid issues with machine precision
        r = suppressWarnings(cor.test(X[, i], X[, j])$estimate)
        if(is.na(r) || r > 0.999)
            next

        if(use_covariates) {
            f_cov = paste(paste(paste0('covariates[, ', 1:ncol(covariates), ']'),
                                collapse=' + '), '+')
        } else {
            f_cov = ''
        }
        f = paste('y ~', f_cov, 'X[, i] * X[, j]')

        glm_result = glm(f, family=fam, na.action=na.omit)
        cf = coefficients(glm_result)[term + 2]
        # Sometimes there is not enough variation to test for
        # an interaction, even after the above filters
        if(is.na(cf))
            next

        p  = coefficients(summary(glm_result))[term + 2, 4]
        N = length(glm_result$fitted.values)
        conv = as.numeric(glm_result$converged)

        cat(chr1, '\t', pos1, '\t', snp1, '\t', maf[i], '\t',
            chr2, '\t', pos2, '\t', snp2, '\t', maf[j], '\t',
            conv, '\t', N,    '\t', p,    '\t', cf,     '\n',
            sep='', file=outhandle)
    }
}
close(outhandle)
